# GETAB_Poldrack.R
# Replication of the investigated effects in a single-subject study
# Date: 2020-10-27
# Author: Alexander Lebedev

# Data from Poldrack et al. (Nature Communications, 2015)
# "Long-term neural and physiological phenotyping of a single human"
# https://www.nature.com/articles/ncomms9885

# Preprocessed data can be downloaded from:
# http://web.stanford.edu/group/poldracklab/myconnectome-data/base/
# OR via the following link (for manuscript review only):
# https://www.dropbox.com/s/pq96x5pa7x5q4zn/PoldrackData.zip?dl=0

# START

# Clear workspace:
rm(list=ls())

# Load libraries:
library(xlsx)
library(tidyverse)
library(quantmod)
library(zoo)
library(xlsx)
library(lmerTest)
library(lme4)
library(nlme)
library(lmtest)
library(reghelper)
library(pracma)
library(tuneR)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(imputeTS)

# Read tracking data:
# Source: http://web.stanford.edu/group/poldracklab/myconnectome-data/base/behavior/trackingdata.txt
poldrackTracking <- read.delim('/Users/alebedev/Documents/Projects/GETAB/data/PoldrackData/trackingdata.txt')
colnames(poldrackTracking)[2] <- 'Date'
poldrackTracking$Date <- as.Date(poldrackTracking$Date)
poldrackTracking[,3:ncol(poldrackTracking)] <- apply(poldrackTracking[,3:ncol(poldrackTracking)],2, as.vector)
poldrackTracking[,3:ncol(poldrackTracking)] <- apply(poldrackTracking[,3:ncol(poldrackTracking)],2, as.numeric)
# SQRT-transform alcohol use patterns:
poldrackTracking$sameevening.AlcoholSquared <- poldrackTracking$sameevening.Alcohol^2
poldrackTracking$prevevening.AlcoholSquared <- poldrackTracking$prevevening.Alcohol^2

# Get market data:
startDate <- range(poldrackTracking$Date)[1]
f = getSymbols("^DJI", auto.assign=FALSE, from=startDate, src='yahoo')
# OR load saved data:
# load('/Users/alebedev/Documents/Projects/GETAB/data/PoldrackData/FTSE100-2020-10-27.rda')
f <- as.data.frame(f)
f$Date <- as.Date(rownames(f))
mName <- names(f)[6] # adjusted values

# Read fMRI data and calculate SD for L and R amygdala:
# Source: http://web.stanford.edu/group/poldracklab/myconnectome-data/base/combined_data_scrubbed/
dat <- merge(poldrackTracking, f, by='Date', all=T)
dat_full <- dat[complete.cases(dat$subcode),] # select complete entries
flist <- list.files('/Users/alebedev/Desktop/Poldr/web.stanford.edu/group/poldracklab/myconnectome-data/base/combined_data_scrubbed/', '*.txt')
fmri <- as.data.frame(matrix(NA, length(flist), 4))
colnames(fmri) <- c('subcode', 'L_Amygdala', 'R_Amygdala', 'LR_AmyCoupling')
for (i in 1:length(flist)){
  fmri$subcode[i] <- strsplit(flist[i],'.txt')[[1]]
  d <- read.delim(paste0('/Users/alebedev/Documents/Projects/GETAB/data/PoldrackData/combined_data_scrubbed/',flist[i]), sep = ' ')
  fmri$L_Amygdala[i] <- 1/sd(d[,622]) # Left Amygdala
  fmri$R_Amygdala[i] <- 1/sd(d[,629]) # Right Amygdala
}

# Merge tracking data with fMRI:
ddd <- merge(dat_full,fmri, by='subcode')

# Remove linear trend:
ddd$DJI.Adjusted.detrended <- detrend(na_interpolation(ddd$DJI.Adjusted))

#################
# Scatter Plots #
#################
# a) systolic blood pressure (proxy for emotional states)
# b) sqrt-transformed alcohol consumption (to achieve close-to-normal distribution)
# c) variance of the BOLD-signal in the Right Amygdala

par(mfrow=c(1,3))
var1= 'afterscan.diastolic'
var2= 'DJI.Adjusted'
cc <- round(cor.test(ddd[,var2], ddd[,var1])$estimate,2)
plot(ddd[,var1], ddd[,var2], pch=16, cex=2,cex.axis=1.5, cex.main=1.5, cex.lab=1.5,
     xlab='Diastolic BP, after scan', ylab='DJI',
     main=paste0('Blood Pressure-Market oscillations: \n Russell Poldrack, r=', cc))
abline(glm(ddd[,var2]~ddd[,var1]), lwd=3)

var1= 'sameevening.AlcoholSquared'
var2= 'DJI.Adjusted'
cc <- round(cor.test(ddd[,var2], ddd[,var1])$estimate,2)
plot(ddd[,var1], ddd[,var2], pch=16, cex=2,cex.axis=1.5, cex.main=1.5, cex.lab=1.5,
     xlab='Alcohol consumption (Squared), same evening', ylab='DJI',
     main=paste0('Alcohol consumption and Market oscillations: \n Russell Poldrack, r=', cc))
abline(glm(ddd[,var2]~ddd[,var1]), lwd=3)

var1= 'DJI.Adjusted'
var2= 'R_Amygdala'
cc <- round(cor.test(ddd[,var2], ddd[,var1])$estimate,2)
plot(ddd[,var1], ddd[,var2], pch=16, cex=4,cex.axis=1.5, cex.main=1.5, cex.lab=1.5,
     xlab='DJI', ylab='1/SD Amygdala (R)',
     main=paste0('Brain-Market oscillations: \n Russell Poldrack, r=', cc))
abline(glm(ddd[,var2]~ddd[,var1]), lwd=3)


par(mfrow=c(1,1))

###########
# TS-plot #
###########
# A few constants
indexColor <- "#69b3a2"
volumeColor <- rgb(0.2, 0.6, 0.9, 1)
dat <- merge(poldrackTracking, f, by='Date', all=T)
ttt <- merge(dat_full,fmri, by='subcode', all=T)
ttt <- ttt[,c(mName, 'Date', 'R_Amygdala')]
# Start from first scan date:
ttt <- ttt[ttt$Date>='2012-10-23',]
# Interpolate Missing Market Data
#ttt[,mName] <- na_interpolation(ttt[,mName])
# Interpolate Missing Scans
#ttt$R_Amygdala <- na_interpolation(ttt$R_Amygdala)

ttt <- ttt[complete.cases(ttt),]
colnames(ttt)[1] <- c('Market')
#ttt$market <-  detrend(tmp_avg$market,tt = 'constant')
coeff=6000
#Plot
ggplot(ttt, aes(x=as.factor(Date))) +
  geom_line( aes(y=Market/2-4000), size=2, color=indexColor, group = 1) + 
  geom_line( aes(y=R_Amygdala*coeff), size=2, color=volumeColor, group = 2) +
  scale_y_continuous(
    name = "Dow Jones Industrial Average",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~., name="Amygdala BOLD-signal variance"),
    #limits=c(3700,8500)
  ) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=7,angle = 90), 
        axis.title.x=element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = indexColor, size=20),
        axis.title.y.right = element_text(color = volumeColor, size=20),
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)
  ) +
  ggtitle("Brain-Market Oscillations: Poldrack data (n=1)")
# END

