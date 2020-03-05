# Create population of representative patients for testing
# from paper Vancomycin PK in normal and morbidly obese subjects
# normal man 1

# ------------------------------------------------------------------------------
# Clear workspace

rm(list=ls(all=TRUE))
graphics.off()

# Set working directory

setwd("D:/Vancomycin")

# Source PopPK model script

source("VancomycinModel.R")

# ------------------------------------------------------------------------------
# Set up an input data frame
# Define values for individual to simulate
set.seed(123456) 
n                  <- 1000                                           # Number of individuals
ID                 <- 1:n                                            # Sequence of individual IDs

AGE                <- 28                                             # Posnatal Age, years
WT                 <- 70.5                                           # Body Weight, kg
SCR                <- 0.914                                          # Serum Creatinine, mg/dL
HEEL               <- 0                                              # Finger-prick Sampling, No = 0
CANCER             <- 0                                              # Cancer, No = 0 
INFDUR             <- 40/60                                          # Infusion Duration, hours

omega.block        <- as.matrix(omat(mod))                           # Omega block values from model
ETA1               <- omega.block[1,1]                               # ETA1 (variance) from model
ETA2               <- omega.block[2,2]                               # ETA2 (variance) from model
ETA3               <- omega.block[3,3]                               # ETA3 (variance) from model

PPVCL              <- rnorm(n,mean = 0,sd = sqrt(ETA1))              # Allocate individuals distribution of PPVCL (SD)
PPVV1              <- rnorm(n,mean = 0,sd = sqrt(ETA2))              # Allocate individuals distribution of PPVV1 (SD)
PPVV2              <- rnorm(n,mean = 0,sd = sqrt(ETA3))              # Allocate individuals distribution of PPVV2 (SD)

patient.data       <- data.frame(ID)
patient.data$AGE   <- AGE
patient.data$WT    <- WT
patient.data$SCR   <- SCR
patient.data$HEEL  <- HEEL
patient.data$CANCER<- CANCER
patient.data$INFDUR<- INFDUR
patient.data$PPVCL <- PPVCL
patient.data$PPVV1 <- PPVV1
patient.data$PPVV2 <- PPVV2
head(patient.data)
