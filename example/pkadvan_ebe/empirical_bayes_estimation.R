# Estimate individual clearance of vemurafenib patients
# -----------------------------------------------------------------------------
# Estimate individual clearance based on concentrations at steady-state. 

# A population pharmacokinetic model for vemurafenib can be found in the FDA
#   Clinical Pharmacology and Biopharmaceutics Review for Vemurafenib.
#   Application Number: 202429Orig1s000 on Page 50. A revised version of this
#   model was proposed by the FDA and can be found on Page 60.

# Ideal method would be to use empirical Bayes estimates from a population
#   pharmacokinetic model. This is not possible using a pharmacokinetic
#   simulation differential equation solver such as `mrgsolve` or `RxODE`.

# Fortunately, as the model is simple we can use analytical solutions
#   instead of differential equations, as the prediction is the 
#   important part, not how we get it.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Source common code
  # source('../../pmg/_source_all.R')

# Read in example data
  source("pkadvan_ebe/pkadvan_simulation.R")

# Load PKADVAN-ready data in NONMEM long format
  # pkdata <- readr::read_rds("datasets/pkadvan_data_tald.rds")
# File had column structure:
# STUDYID - study identifier (used with ddply, but just SUBJID would suffice)
# SUBJID - subject identifier (used with ddply)
# VISIT - visit identifier (not used by EBE function)
# TIME - time after first dose
# AMT - dose
# MDV - missing dependent variable
# DV - dependent variable
# TALD - time after last dose (not used by EBE function)
# DAY - days starting on day of first dose
# STDY - study day (not used by EBE function)
# FDTIME - time of first dose in hms format (not used by EBE function)
# FDDAY - study day of first dose (not used by EBE function)
# WT - weight

# Load required libraries
  library(PKADVAN)  # ??PKADVAN for vignette on how PKADVAN works

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
# Define empirical bayesian estimation process function
  ebe_fn <- function(input) {
  # Prepare environment
  # Convert input to data.table format
    input <- as.data.table(input)[TIME >= 0, ]
  # Define PPV vector and inital values for ETAs
    ETAPPV <- with(pk_mod, c(PPVCL, PPVV, PPVKA))  # CV is the desired format
    n_eta <- length(ETAPPV)
    init_par <- exp(double(n_eta))
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Set up loop to repeat bayesian estimation until predictions for individual
  #   parameters minimise successfully.
    run_num <- 0
    repeat {
    # If there has been any previous unsuccessful runs, sample initial ETAs
    #   from random uniform distribution
      if (run_num > 0) {
        init_par <- init_par*exp(runif(n_eta, min = -0.01, max = 0.01))
      }
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Define bayesian estimation function
      bayes_estimate_fn <- function(par) {
      # Describe parameters to be optimised and dataset to be used
        ETA <- log(par)
        pk_in <- input
      # Define PKADVAN input dataset
      # Relative bioavailability - time-varying covariate
        pk_in[DAY <= 14 & STUDYID == "BRIM2", F1 := with(pk_mod, POPF1*Tlt14onF1)]
        pk_in[DAY > 14 & DAY <= 104 & STUDYID == "BRIM2", F1 := with(pk_mod, POPF1*Tlt104onF1)]
        pk_in[(DAY >= 105 & STUDYID == "BRIM2") | STUDYID != "BRIM2", F1 := with(pk_mod, POPF1)]
      # Clearance and volume - covariate and population parameter variability
        pk_in[, CL := with(pk_mod, POPCL * (WT/REFWT)^WTonCL * exp(ETA[1]) )]
        pk_in[, V := with(pk_mod, POPV * (WT/REFWT)^WTonV * exp(ETA[2]) )]
      # Absorption constant - population parameter variability
        pk_in[, KA := with(pk_mod, POPKA * exp(ETA[3]) )]
      # Calculate concentrations using PKADVAN
        pk_out <- PKADVAN::OneCompFirstOrderAbs(pk_in)
      # Ensure DV has finite values
        pk_out[!is.finite(IPRED) | IPRED < .Machine$double.eps, IPRED := .Machine$double.eps]
      # Define DV (y) and individual prediction for DV (yhat)
        y <- pk_out[MDV == 0, DV]
        yhat <- pk_out[MDV == 0, IPRED]
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Determine posterior and prior log-likelihood
      # Posterior log-likelihood
      # Error model: IPRED*(1+ + RUVPRO*EPS(1)) + RUVADD*(EPS(2))
      # Can be simplified to: IPRED + W*EPS(1)
      # where W = sqrt((IPRED*RUVPRO)^2 + RUVADD^2)
        loglikpost_sd <- with(pk_mod, sqrt((yhat*RUVPRO)^2 + RUVADD^2))
        loglikpost <- dnorm(y, mean = yhat, sd = loglikpost_sd, log = T)
      # Prior log-likelihood
        loglikprior <- dnorm(ETA, mean = 0, sd = ETAPPV, log = T)
      # Return objective function value to be minimised
        return(-1*sum(loglikpost, loglikprior))
      }
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Use optim to optimise the objective function returned by bayes_estimate_fn
      bayes_estimate <- try(optim(init_par, bayes_estimate_fn, method = "L-BFGS-B",
        lower = rep(0.001, times = n_eta), upper = rep(1000, times = n_eta),
        control = list(
          parscale = init_par, fnscale = bayes_estimate_fn(init_par),
          factr = 1e12, pgtol = 1e-8
        )
      ))
    # Check to see if the minimisation was successful, if not repeat
      if (class(bayes_estimate) == "try-error") browser()  # error catching
      minimised <- "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
      if (bayes_estimate$message == minimised) break  # exits repeat
      run_num <- run_num + 1
      if (run_num > 10) browser()  # multiple unsuccessful runs
    }  # end repeat
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Use empirical bayes estimates to give desired output
  # While we mainly want the empirical bayes estimates for clearance, we will
  #   also collect IPRED and PRED to assess the model's ability to model data
  #   from coBRIM which wasn't used for model development.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Individual Predicted
    pkadvan_ipred <- input
  # Extract final ETA values from our minimised bayesian estimation step
    ETA <- log(bayes_estimate$par)
    ETA_tbl <- as.data.table(matrix(rep(ETA, dim(input)[1]), ncol = n_eta, byrow = TRUE))
    names(ETA_tbl) <- paste0("ETA", readr::parse_number(names(ETA_tbl)))
  # Define PKADVAN input dataset
  # Relative bioavailability - time-varying covariate
    pkadvan_ipred[DAY <= 14 & STUDYID == "BRIM2", F1 := with(pk_mod, POPF1*Tlt14onF1)]
    pkadvan_ipred[DAY > 14 & DAY <= 104 & STUDYID == "BRIM2", F1 := with(pk_mod, POPF1*Tlt104onF1)]
    pkadvan_ipred[(DAY >= 105 & STUDYID == "BRIM2") | STUDYID != "BRIM2", F1 := with(pk_mod, POPF1)]
  # Clearance and volume - covariate and population parameter variability
    pkadvan_ipred[, CL := with(pk_mod, POPCL * (WT/REFWT)^WTonCL * exp(ETA[1]) )]
    pkadvan_ipred[, V := with(pk_mod, POPV * (WT/REFWT)^WTonV * exp(ETA[2]) )]
  # Absorption constant - population parameter variability
    pkadvan_ipred[, KA := with(pk_mod, POPKA * exp(ETA[3]) )]
  # Predict individual concentrations using PKADVAN
    ipred_out <- PKADVAN::OneCompFirstOrderAbs(pkadvan_ipred)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Population Predicted
  # We change the ETA's to zero (i.e. log(1)) and leave the rest of the model
  #   the same.
    pkadvan_pred <- input
  # Extract final ETA values from our minimised bayesian estimation step
    ETA <- log(rep(1, n_eta))  # shrink ETA to zero
  # Define PKADVAN input dataset (same code)
  # Relative bioavailability - time-varying covariate
    pkadvan_pred[DAY <= 14 & STUDYID == "BRIM2", F1 := with(pk_mod, POPF1*Tlt14onF1)]
    pkadvan_pred[DAY > 14 & DAY <= 104 & STUDYID == "BRIM2", F1 := with(pk_mod, POPF1*Tlt104onF1)]
    pkadvan_pred[(DAY >= 105 & STUDYID == "BRIM2") | STUDYID != "BRIM2", F1 := with(pk_mod, POPF1)]
  # Clearance and volume - covariate and population parameter variability
    pkadvan_pred[, CL := with(pk_mod, POPCL * (WT/REFWT)^WTonCL * exp(ETA[1]) )]
    pkadvan_pred[, V := with(pk_mod, POPV * (WT/REFWT)^WTonV * exp(ETA[2]) )]
  # Absorption constant - population parameter variability
    pkadvan_pred[, KA := with(pk_mod, POPKA * exp(ETA[3]) )]
  # Predict individual concentrations using PKADVAN
    pred_out <- PKADVAN::OneCompFirstOrderAbs(pkadvan_pred)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return desired output
    output <- cbind(ipred_out, ETA_tbl)
    output$PRED <- pred_out$IPRED  # pred_out$IPRED is actually PRED
  # Bind in true values for ETAs  
  # Calculate weighted residuals and individual weighted residuals
    output[, WRES := (DV-PRED)/with(pk_mod, sqrt((PRED*RUVPRO)^2 + RUVADD^2))]
    output[, IWRES := (DV-PRED)/with(pk_mod, sqrt((IPRED*RUVPRO)^2 + RUVADD^2))]
  # Cal
    return(output)
  }
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Estimate empirical bayes estimates for each subject
# WARNING: This takes a significant amount of time, you might want to grab a 
#   cup of tea, watch a long youtube video, add comments to your code, do
#   some other work or something.
  ebedata <- as.data.table(ddply(pkdata, .(STUDYID, SUBJID), ebe_fn))
  
# Save EBE data for later use
  readr::write_rds(ebedata, "pkadvan_simulation/ebe_data.rds")
  