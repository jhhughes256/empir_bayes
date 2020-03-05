# Busulfan - Create population of representative patients for regimen testing
# Westmead Example Patient

# ------------------------------------------------------------------------------
# Clear workspace
rm(list=ls(all=TRUE))
graphics.off()

# Set working directory
setwd("D:/Busulfan")

# # Load pakage libraries
library(dplyr)	    #New plyr - required for mrgsolve
library(mrgsolve)	  #Metrum differential equation solver for pharmacometrics
# 
# Source PopPK model script
source("BusulfanModel.R")

# ------------------------------------------------------------------------------
# Set up an input data frame
# Define values for individual to simulate
set.seed(123456) 
n           <- 1000                                                              # Number of individuals
ID          <- 1:n                                                               # Sequence of individual ID's

BW          <- 71.2                                                                # Weight (kg)
SEX         <- 0                                                                 # Sex (0 = male, 1 = female)
HEIGHT      <- 150                                                               # Height (cm)
INFDUR      <- 200/60                                                                # Infusion Duration, hours


omega.block <- as.matrix(omat(mod))                 # Omega block values from model
ETA1       <- omega.block[1,1]                     # PPVCL (variance) from model
ETA2       <- omega.block[2,2]                     # PPVV2 (variance) from model


PPVCL        <- rnorm(n,mean = 0,sd = sqrt(ETA1))   # Allocate individuals distribution of PPVCL (SD)
PPVV1        <- rnorm(n,mean = 0,sd = sqrt(ETA2))   # Allocate individuals distribution of PPVV1 (SD)


# Create data frame of individuals with varying demographics and ETA values
patient.data                                     <- data.frame(ID)
patient.data$BW                                  <- BW
patient.data$SEX                                 <- SEX
patient.data$HEIGHT                              <- HEIGHT
patient.data$INFDUR                              <- INFDUR

patient.data$PPVCL                               <- PPVCL
patient.data$PPVV1                               <- PPVV1

head(patient.data)

