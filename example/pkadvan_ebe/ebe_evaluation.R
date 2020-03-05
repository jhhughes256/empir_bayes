# Evaluation of empirical bayesian estimates
# -----------------------------------------------------------------------------
# Evaluate the empirical bayesian estimates provided by the script:
#  - empirical_bayes_estimation.R

# Estimates come from vemurafenib model:
# pk_mod <- list(
#   # Population Parameters
#   POPCL = 1.3,  # apparent clearance; L/h
#   POPV = 106,  # apparent volume; L
#   POPKA = 0.188,  # absorption constant; h-1
#   POPF1 = 1,  # relative bioavailability, DAY >= 105
#   # Covariate Effects
#   REFWT = 70,
#   WTonCL = 0.319,
#   WTonV = 0.740,
#   Tlt14onF1 = 0.789,  # relative bioavailability, DAY <= 14
#   Tlt104onF1 = 0.899,  # relative bioavailability, DAY > 15 & DAY <= 104
#   # Population Parameter Variability
#   PPVCL = 0.319,  # apparent clearance variability; CV (coefficient of variation)
#   PPVV = 0.657,  # apparent volume variability; CV
#   PPVKA = 1.01,  # absorption constant variability; CV
#   # Residual Unexplained Variability
#   RUVADD = 0.814,  # additive error; mg/L
#   RUVPRO = 0.228  # proportional error; CV
# )

# Evaluation includes:
# - Individual parameter distribution
# - ETA Shrinkage
# - Individual Predicted vs. Observed
# - Weighted residual plots
# - VPC

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Source common code
  # source('../../pmg/_source_all.R')

# Load PKADVAN-ready data in NONMEM long format
  source("pkadvan_ebe/empirical_bayes_estimation.R")
  # ebedata <- readr::read_rds("pkadvan_ebe/ebe_data.rds")
  
# Set colourblind palette
  cbPalette <- data.frame(
    grey = "#999999",
    orange = "#E69F00",
    skyblue = "#56B4E9",
    green = "#009E73",
    yellow = "#F0E442",
    blue = "#0072B2",
    red = "#D55E00",
    pink = "#CC79A7",
    stringsAsFactors = F
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Original model was built with BRIM2, BRIM3 and unavailable Phase I data
# Determine whether coBRIM is modelled well when assumed to be like BRIM3
  ebedata[!duplicated(SUBJID) & STUDYID == "coBRIM", describe(CL)]
  ebedata[!duplicated(SUBJID) & STUDYID == "coBRIM", plot(density(log(CL)))]
  
# Overlay estimated clearance values over distribution
  nid <- ebedata[STUDYID == "coBRIM", length(unique(SUBJID))] - 2
  ind_wt <- ebedata[!duplicated(SUBJID) & STUDYID == "coBRIM", WT]
  cl_sim_cobrim <- 1.3 * (rep(ind_wt, each = 1000)/70)^0.319 * exp(rnorm(1000*nid, 0, 0.319))
  plot(density(log(cl_sim_cobrim)))

  p <- NULL
  p <- ggplot() + theme_bw()
  p <- p + geom_boxplot(aes(x = "SIM", y = cl_sim_cobrim), outlier.shape = NA)
  p <- p + geom_boxplot(aes(x = "EBE", 
    y = ebedata[!duplicated(SUBJID) & STUDYID == "coBRIM", CL]), outlier.shape = NA)
  p <- p + geom_dotplot(aes(x = "EBE", 
    y = ebedata[!duplicated(SUBJID) & STUDYID == "coBRIM", CL]),
    stackdir = "center", binaxis = "y", binwidth = 0.05, alpha = 0.1)
  p <- p + labs(x = NULL, y = "Estimated Clearance (L/h)")
  p_boxdot <- p + coord_cartesian(xlim = NULL, ylim = c(0.3, 4.4))
  p_boxdot
  
  p <- NULL
  p <- ggplot() + theme_bw()
  p <- p + geom_density(aes(x = cl_sim_cobrim), size = 1, colour = cbPalette$blue)
  p <- p + geom_vline(xintercept = median(cl_sim_cobrim), alpha = 0.5, size = 1, colour = cbPalette$blue)
  p <- p + geom_vline(xintercept = quantile(cl_sim_cobrim, probs = c(0.05, 0.95)), alpha = 0.5, 
    size = 1, linetype = "dashed", colour = cbPalette$blue)
  p <- p + geom_density(aes(x = ebedata[!duplicated(SUBJID) & STUDYID == "coBRIM", CL]),
    colour = cbPalette$red, size = 1)
  p <- p + geom_vline(xintercept = rep(median(ebedata[!duplicated(SUBJID) & STUDYID == "coBRIM", CL]), 2),
    alpha = 0.5, size = 1, colour = cbPalette$red)
  p <- p + geom_vline(xintercept = quantile(ebedata[!duplicated(SUBJID) & STUDYID == "coBRIM", CL], probs = c(0.05, 0.95)),
    alpha = 0.5, size = 1, linetype = "dashed", colour = cbPalette$red)
  p <- p + labs(x = "Clearance (L/h)", y = "Density")
  p_dens <- p + coord_cartesian(xlim = NULL, ylim = c(0.055, 1.15))
  p_dens
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Now do the same but for all data, not just coBRIM
  ebedata[!duplicated(SUBJID), describe(CL)]
  ebedata[!duplicated(SUBJID), plot(density(log(CL)))]
  
# Overlay estimated clearance values over distriution
  nid <- ebedata[, length(unique(SUBJID))]
  ind_wt <- ebedata[!duplicated(SUBJID), WT]
  cl_sim <- 1.3 * (rep(ind_wt, each = 1000)/70)^0.319 * exp(rnorm(1000*nid, 0, 0.319))
  plot(density(log(cl_sim)))
  
  p <- NULL
  p <- ggplot() + theme_bw()
  p <- p + geom_boxplot(aes(x = "SIM", y = cl_sim), outlier.shape = NA)
  p <- p + geom_boxplot(aes(x = "EBE", 
    y = ebedata[!duplicated(SUBJID), CL]), outlier.shape = NA)
  p <- p + geom_dotplot(aes(x = "EBE", 
    y = ebedata[!duplicated(SUBJID), CL]),
    stackdir = "center", binaxis = "y", binwidth = 0.04, alpha = 0.1)
  p <- p + labs(x = NULL, y = "Estimated Clearance (L/h)")
  p_boxdot <- p + coord_cartesian(xlim = NULL, ylim = c(0.3, 4.4))
  p_boxdot
  
  p <- NULL
  p <- ggplot() + theme_bw()
  p <- p + geom_density(aes(x = cl_sim), size = 1, colour = cbPalette$blue)
  p <- p + geom_vline(xintercept = median(cl_sim), alpha = 0.5, size = 1, colour = cbPalette$blue)
  p <- p + geom_vline(xintercept = quantile(cl_sim, probs = c(0.05, 0.95)), 
    alpha = 0.5, size = 1, linetype = "dashed", colour = cbPalette$blue)
  p <- p + geom_density(aes(x = ebedata[!duplicated(SUBJID), CL]),
    colour = cbPalette$red, alpha = 0.6, size = 1)
  p <- p + geom_vline(xintercept = rep(median(ebedata[!duplicated(SUBJID), CL]), 2),
    alpha = 0.5, size = 1, colour = cbPalette$red)
  p <- p + geom_vline(xintercept = quantile(ebedata[!duplicated(SUBJID), CL], probs = c(0.05, 0.95)),
    alpha = 0.5, size = 1, linetype = "dashed", colour = cbPalette$red)
  p <- p + labs(x = "Clearance (L/h)", y = "Density")
  p_dens <- p + coord_cartesian(xlim = NULL, ylim = c(0.055, 1.15))
  p_dens
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ETA Shrinkage
# Shrinkage implies bias, ideally it is below 20-30%.
# High ETA shrinkage implies that predictions are biased to shrink towards the
#   mean value. High ETA shrinkage will result in EBEs that are not indicative
#   of the individual, but instead the population.
# ETA shrinkage: 1 - SD(ETA)/PPV
  list(
    CL_ETA_shrinkage = 1 - sd(ebedata[!duplicated(SUBJID), ETA1])/0.319,
    Vd_ETA_shrinkage = 1 - sd(ebedata[!duplicated(SUBJID), ETA2])/0.657,
    ka_ETA_shrinkage = 1 - sd(ebedata[!duplicated(SUBJID), ETA3])/1.01
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Individual vs Observed - External Model Validation
  pdata <- ebedata[MDV == 0, ]
  pdata[TALD < 0, TALD := TALD + 12]
  myPalette <- with(cbPalette, c(red, blue, green))
  
  p <- NULL
  p <- ggplot(aes(x = PRED, y = DV), data = pdata)
  p <- p + theme_bw()
  p <- p + geom_point(aes(colour = STUDYID), size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + scale_colour_manual("Study", values = myPalette)
  p <- p + labs(x = "Predicted (mg/L)", y = "Observed (mg/L)")
  p_DVvPRED <- p
  p_DVvPRED
  
  p <- NULL
  p <- ggplot(aes(x = IPRED, y = DV), data = pdata)
  p <- p + theme_bw()
  p <- p + geom_point(aes(colour = STUDYID), size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + scale_colour_manual("Study", values = myPalette)
  p <- p + labs(x = "Individual Predicted (mg/L)", y = "Observed (mg/L)")
  p_DVvIPRE <- p
  p_DVvIPRE
  
  p <- NULL
  p <- ggplot(aes(x = TIME/24, y = WRES), data = pdata)
  p <- p + theme_bw()
  p <- p + geom_point(aes(colour = STUDYID), size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + scale_colour_manual("Study", values = myPalette)
  p <- p + labs(x = "Weighted Residual", y = "Time (days)")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(-12.5, 12.5))
  p_TIMEvWRES <- p
  p_TIMEvWRES
  
  p <- NULL
  p <- ggplot(aes(x = TALD, y = WRES), data = pdata)
  p <- p + theme_bw()
  p <- p + geom_point(aes(colour = STUDYID), size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + scale_colour_manual("Study", values = myPalette)
  p <- p + labs(x = "Time after last dose (hours)", y = "Weighted Residual")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(-12.5, 12.5))
  p_TALDvWRES <- p
  p_TALDvWRES
  
  p <- NULL
  p <- ggplot(aes(x = PRED, y = WRES), data = pdata)
  p <- p + theme_bw()
  p <- p + geom_point(aes(colour = STUDYID), size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + scale_colour_manual("Study", values = myPalette)
  p <- p + labs(x = "Predicted (mg/L)", y = "Weighted Residual")
  p <- p + coord_cartesian(xlim = c(0, 100), ylim = c(-12.5, 12.5))
  p_PREDvWRES <- p
  p_PREDvWRES
  
# Check coBRIM specifically
  p <- NULL
  p <- ggplot(aes(x = PRED, y = DV), data = pdata[STUDYID == "coBRIM", ])
  p <- p + theme_bw()
  p <- p + geom_point(colour = cbPalette$green, size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + labs(x = "Predicted (mg/L)", y = "Observed (mg/L)")
  p_DVvPRED_coBRIM <- p
  p_DVvPRED_coBRIM
  
  p <- NULL
  p <- ggplot(aes(x = IPRED, y = DV), data = pdata[STUDYID == "coBRIM", ])
  p <- p + theme_bw()
  p <- p + geom_point(colour = cbPalette$green, size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + labs(x = "Individual Predicted (mg/L)", y = "Observed (mg/L)")
  p_DVvIPRE_coBRIM <- p
  p_DVvIPRE_coBRIM
  
  p <- NULL
  p <- ggplot(aes(x = TIME/24, y = WRES), data = pdata[STUDYID == "coBRIM", ])
  p <- p + theme_bw()
  p <- p + geom_point(colour = cbPalette$green, size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + scale_colour_manual("Study", values = myPalette)
  p <- p + labs(x = "Time (days)", y = "Weighted Residual")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(-12.5, 12.5))
  p_TIMEvWRES_coBRIM <- p
  p_TIMEvWRES_coBRIM
  
  p <- NULL
  p <- ggplot(aes(x = TALD, y = WRES), data = pdata[STUDYID == "coBRIM", ])
  p <- p + theme_bw()
  p <- p + geom_point(colour = cbPalette$green, size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + scale_colour_manual("Study", values = myPalette)
  p <- p + labs(x = "Time after last dose (hours)", y = "Weighted Residual")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(-12.5, 12.5))
  p_TALDvWRES_coBRIM <- p
  p_TALDvWRES_coBRIM
  
  p <- NULL
  p <- ggplot(aes(x = PRED, y = WRES), data = pdata[STUDYID == "coBRIM", ])
  p <- p + theme_bw()
  p <- p + geom_point(colour = cbPalette$green, size = 1.5, shape = 1, alpha = 0.5)
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + stat_smooth(colour = "red", size = 1, method = "loess")
  p <- p + scale_colour_manual("Study", values = myPalette)
  p <- p + labs(x = "Predicted (mg/L)", y = "Weighted Residual")
  p <- p + coord_cartesian(xlim = NULL, ylim = c(-12.5, 12.5))
  p_PREDvWRES_coBRIM <- p
  p_PREDvWRES_coBRIM
  
# Could also do VPC, but thats a bit overkill I think