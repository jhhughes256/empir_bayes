# Simulate PKADVAN data for empirical bayes estimation example
# -----------------------------------------------------------------------------
# Prepare work environment
# Load required libraries
  library(PKADVAN)  # ??PKADVAN for vignette on how PKADVAN works
  library(data.table)
  library(plyr)

# Define pharmacokinetic model
# The pharmacokinetics of vemurafenib are explained by a one-compartment model
#   with first order absorption & first order elimination.
# Model was built using BRIM2, BRIM3 and an unavailable Phase I dataset
# PKADVAN function for this is: PKADVAN::OneCompFirstOrderAbs
# Requires parameters for F1, KA, CL & V (found using ?PKADVAN::OneCompFirstOrderAbs)
# Order and naming doesn't matter *here*, just a way of keeping our model 
#   parameters tidy and not scattered into the global environment.
  pk_mod <- list(
  # Population Parameters
    POPCL = 1.3,  # apparent clearance; L/h
    POPV = 106,  # apparent volume; L
    POPKA = 0.188,  # absorption constant; h-1
    POPF1 = 1,  # relative bioavailability, DAY >= 105
  # Covariate Effects
    REFWT = 70,
    WTonCL = 0.319,
    WTonV = 0.740,
    Tlt14onF1 = 0.789,  # relative bioavailability, DAY <= 14  # only on BRIM2
    Tlt104onF1 = 0.899,  # relative bioavailability, DAY > 15 & DAY <= 104  # only on BRIM2
  # Population Parameter Variability
    PPVCL = 0.319,  # apparent clearance variability; CV (coefficient of variation)
    PPVV = 0.657,  # apparent volume variability; CV
    PPVKA = 1.01,  # absorption constant variability; CV
  # Residual Unexplained Variability
    RUVADD = 0.814,  # additive error; mg/L
    RUVPRO = 0.228  # proportional error; CV
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create subject level data
# Define subject numbers
  nid <- 60
  ID <- 1:nid
  
# Sample Individual Population Parameters
  pop_tbl <- data.table(
    STUDYID = rep(c("BRIM2", "BRIM3", "coBRIM"), each = 20),
    SUBJID = ID,
    ETA1 = rnorm(nid, 0, sd = sqrt(pk_mod$PPVCL)),
    ETA2 = rnorm(nid, 0, sd = sqrt(pk_mod$PPVV)),
    ETA3 = rnorm(nid, 0, sd = sqrt(pk_mod$PPVKA)),
    WT = rnorm(nid, 70, sd = 5)
  )
  
# Create longitudinal data
# Dose data
  dose_times <- seq(0, 1020, by = 12)
  dose_tbl <- data.table(
    SUBJID = rep(ID, each = length(dose_times)),
    TIME = dose_times,
    AMT = 960,
    MDV = 1,
    DV = 0
  )
    
# Concentration data
  conc_times_brim2 <- as.numeric(outer(c(-0.25, 2, 4, 6, 8), c(0, 360, 504, 1008), FUN = "+"))
  conc_brim2 <- data.table(
    SUBJID = rep(ID[1:20], each = length(conc_times_brim2)),
    TIME = rep(conc_times_brim2, nid/3),
    AMT = 0,
    MDV = 0,
    DV = 0
  )
  conc_times_brim3 <- as.numeric(outer(c(-0.25, 2.5), c(0, 360, 504), FUN = "+"))
  conc_brim3 <- data.table(
    SUBJID = rep(ID[21:40], each = length(conc_times_brim3)),
    TIME = rep(conc_times_brim3, nid/3),
    AMT = 0,
    MDV = 0,
    DV = 0
  )
  conc_times_cobrim <- as.numeric(outer(c(-0.25, 3), c(0, 360, 672), FUN = "+"))
  conc_cobrim <- data.table(
    SUBJID = rep(ID[41:60], each = length(conc_times_cobrim)),
    TIME = rep(conc_times_cobrim, nid/3),
    AMT = 0,
    MDV = 0,
    DV = 0
  )
  conc_tbl <- rbind(conc_brim2, conc_brim3, conc_cobrim)
  conc_tbl <- conc_tbl[TIME >= 0, ]
  
# Merge final dataset
  pk_tbl <- rbind(dose_tbl, conc_tbl)
  input_tbl <- merge(pop_tbl, pk_tbl, by = "SUBJID")
  setorder(input_tbl, SUBJID, TIME)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define simulation function
  pkadvan_sim_fn <- function(input_df) {
  # Define parameter values
    ipred_in <- as.data.table(input_df)
    ipred_in[, DAY := TIME/24 + 1]
    ipred_in[, F1 := with(pk_mod, POPF1)]
    ipred_in[DAY <= 14 & STUDYID == "BRIM2", F1 := with(pk_mod, POPF1*Tlt14onF1)]
    ipred_in[DAY > 14 & STUDYID == "BRIM2", F1 := with(pk_mod, POPF1*Tlt104onF1)]
  # Clearance and volume - covariate and population parameter variability
    ipred_in[, CL := with(pk_mod, POPCL * (WT/REFWT)^WTonCL * exp(ETA1) )]
    ipred_in[, V := with(pk_mod, POPV * (WT/REFWT)^WTonV * exp(ETA2) )]
  # Absorption constant - population parameter variability
    ipred_in[, KA := with(pk_mod, POPKA * exp(ETA3) )]
  # Predict individual concentrations using PKADVAN
    ipred_out <- PKADVAN::OneCompFirstOrderAbs(ipred_in)
  # Sample residual unexplained variability
    ipred_out[, W := with(pk_mod, sqrt((IPRED*RUVPRO)^2) + RUVADD^2)]
    ipred_out[, DV := IPRED + W*rnorm(length(IPRED), 0, 1)]
  # Clean up and return data
    ipred_out[DV <= 0, DV := 0]
    ipred_out[DV <= 0, MDV := 1]
    ipred_out[, c("STUDYID", "SUBJID", "TIME", "AMT", "MDV", "DV", "WT", "DAY")]
  }
  pkdata <- as.data.table(ddply(input_tbl, .(SUBJID), pkadvan_sim_fn))
  