# Simulate representative patients with dosing regimen
# first normal man in paper "Vancomycin PK in normal and morbidly obese subjects"

# ------------------------------------------------------------------------------
# Clear workspace

rm(list=ls(all=TRUE))
graphics.off()

# Set working directory

setwd("D:/Vancomycin")

# Source PopPK model script

source("VancomycinPopulation.R")

# ------------------------------------------------------------------------------
# Replicate test population for concentration dataset
TIME.conc          <- seq(from = 0, to = 24, by = 5/60)
conc.data          <- lapply(patient.data,rep.int, times=length(TIME.conc))
conc.data          <- as.data.frame(conc.data)
conc.data          <- conc.data[with(conc.data, order(conc.data$ID)),]
conc.data$TIME     <- TIME.conc
conc.data$AMT      <- 0
conc.data$CMT      <- 1
conc.data$EVID     <- 0
conc.data$RATE     <- 0
head(conc.data)

# Replicate test population for dose dataset
TIME.dose          <- c(0)
dose.data          <- lapply(patient.data,rep.int, times=length(TIME.dose))
dose.data          <- as.data.frame(dose.data)
dose.data          <- dose.data[with(dose.data, order(dose.data$ID)),]
dose.data$TIME     <- TIME.dose
dose.data$AMT      <- 1000       # dose in mg
dose.data$CMT      <- 1
dose.data$EVID     <- 1
dose.data$RATE     <- -2        
head(dose.data)

# Combine into simulation dataset
input.data         <- rbind(conc.data, dose.data)
input.data         <- input.data[with(input.data, order(input.data$ID, input.data$TIME, input.data$EVID)),]
head(input.data)

# ------------------------------------------------------------------------------
# Simulate dataset
output.data        <- mod %>% data_set(input.data) %>% carry.out(AMT, EVID) %>% mrgsim()
output.data        <- as.data.frame(output.data)
head(output.data)

#Save simulated dataset
# write.csv(output.data,file="VancomycinPopPKSimulationData.csv")

#------------------------------------------------------------------------------------
# Plot simulation results
# Use a custom ggplot2 theme
theme_bw2          <- theme_set(theme_bw(base_size = 12))
theme_bw2          <- theme_update(plot.margin = unit(c(1,1,1,1), "lines"),
                                   axis.title.x=element_text(size = 12, vjust = 0),
                                   axis.title.y=element_text(size = 12, vjust = 0, angle = 90),
                                   strip.text.x=element_text(size = 10),
                                   strip.text.y=element_text(size = 10, angle = 90))

#Define 90% CI functions
CI90lo             <- function(x) quantile(x, probs = 0.05)
CI90hi             <- function(x) quantile(x, probs = 0.95)

#Subset data to remove dosing events
plot.data          <- subset(output.data,EVID==0)
head(plot.data)

#Create concentration-time plot
conctime.plot      <- NULL
conctime.plot      <- ggplot()
conctime.plot      <- conctime.plot + stat_summary(aes(x=TIME,y=IPRE),data=plot.data,geom="line",fun.y=median,colour="dodgerblue4",size=1)
conctime.plot      <- conctime.plot + stat_summary(aes(x=TIME,y=IPRE),data=plot.data,geom="ribbon",fun.ymin="CI90lo",fun.ymax="CI90hi",fill="dodgerblue4",alpha=0.2)
conctime.plot      <- conctime.plot + scale_y_continuous("Vancomycin Concentration (mg/dL)\n")
conctime.plot      <- conctime.plot + theme(axis.title.y = element_text(face="bold", angle=90, size=14),axis.text.y=element_text(size=10))
conctime.plot      <- conctime.plot + scale_x_continuous("\nTime (hr)")
conctime.plot      <- conctime.plot + theme(axis.title.x = element_text(face="bold", size=14),axis.text.x=element_text(size=10))
print(conctime.plot)

#Save plot
# ggsave("VancomycinPopPKSimulationPlot.png", height = 8, width = 8)