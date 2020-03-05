# Busulfan  -  Simulate representative patient population with test dosing regimen
# Westmead Example Patient
# 3.2mg/kg iv

# ------------------------------------------------------------------------------
# Clear workspace
rm(list=ls(all=TRUE))
graphics.off()

# Set working directory
setwd("D:/Busulfan")

# Load pakage libraries
library(plyr)	      #New plyr - required for mrgsolve
library(dplyr)	    #Split and rearrange data - required for mrgsolve
library(mrgsolve)	  #Metrum differential equation solver for pharmacometrics
library(ggplot2)    #Graphical package

# Source PopPK model script
source("BusulfanPopulation.R")    #with infusion duration of 3.33 hours


# ------------------------------------------------------------------------------
# Replicate test population for concentration dataset
TIME.conc           <- seq(from = 0, to = 24, by = 0.33)
test.conc           <- lapply(patient.data,rep.int, times=length(TIME.conc))
test.conc           <- as.data.frame(test.conc)
test.conc           <- test.conc[with(test.conc, order(test.conc$ID)),]
test.conc$time      <- TIME.conc
test.conc$amt       <- 0
test.conc$cmt       <- 1
test.conc$evid      <- 0
test.conc$rate      <- 0 #duration code
head(test.conc)

# Replicate test population for Day 1 dose dataset
TIME.dose1          <- c(0)
test.dose1          <- lapply(patient.data,rep.int, times=length(TIME.dose1))
test.dose1          <- as.data.frame(test.dose1)
test.dose1          <- test.dose1[with(test.dose1, order(test.dose1$ID)),]
test.dose1$time     <- TIME.dose1
test.dose1$amt      <- 3.2*test.dose1$BW
test.dose1$cmt      <- 1
test.dose1$evid     <- 1
test.dose1$rate     <- -2 #duration code
head(test.dose1)

# Combine into simulation dataset
test.input          <- rbind(test.conc, test.dose1)
test.input          <- test.input[with(test.input, order(test.input$ID, test.input$time, test.input$evid)),]
head(test.input)

# ------------------------------------------------------------------------------
# Simulate standard dosing regimen
test.data           <- mod %>% data_set(test.input) %>% carry.out(amt, evid) %>% mrgsim()
test.data           <- as.data.frame(test.data)
head(test.data)

#Save simulated dataset
# write.csv(test.data,file="BusulfanPopPKSimulationData.csv")

#------------------------------------------------------------------------------------
#Plot simulation results
#Use a custom ggplot2 theme
theme_bw2           <- theme_set(theme_bw(base_size = 12))
theme_bw2           <- theme_update(plot.margin = unit(c(1,1,1,1), "lines"),
                                    axis.title.x=element_text(size = 12, vjust = 0),
                                    axis.title.y=element_text(size = 12, vjust = 0, angle = 90),
                                    strip.text.x=element_text(size = 10),
                                    strip.text.y=element_text(size = 10, angle = 90))

#Define 90% CI functions
CI90lo              <- function(x) quantile(x, probs = 0.05)
CI90hi              <- function(x) quantile(x, probs = 0.95)

#Subset data to remove dosing events
test.data2          <- subset(test.data,evid==0)
AUC.data            <- subset(test.data,time==24)

#Create concentration-time plot
test.plot           <- NULL	
test.plot           <- ggplot()	
test.plot           <- test.plot + stat_summary(aes(x=time,y=IPRE),data=test.data2,geom="line",fun.y=median,colour="indianred3",size=1)	
test.plot           <- test.plot + stat_summary(aes(x=time,y=IPRE),data=test.data2,geom="ribbon",fun.ymin="CI90lo",fun.ymax="CI90hi",fill="indianred3",alpha=0.2)
test.plot           <- test.plot + scale_y_continuous("Busulfan Concentration (mg/L)\n")
test.plot           <- test.plot + theme(axis.title.y = element_text(face="bold", angle=90, size=18),axis.text.y=element_text(size=16))
test.plot           <- test.plot + scale_x_continuous("\nTime (hr)")
test.plot           <- test.plot + theme(axis.title.x = element_text(face="bold", size=18),axis.text.x=element_text(size=16))
print(test.plot)	

#Save plot
# ggsave("BusulfanPopPKSimulationPlot.png", height = 8, width = 8) 