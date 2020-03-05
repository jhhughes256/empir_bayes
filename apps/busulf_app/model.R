# Salinger 2010, Busulfan Population Pharmacokinetic Model
# ------------------------------------------------------------------------------
# Load libraries
  # library(dplyr)	    #New plyr - required for mrgsolve
  # library(mrgsolve)	  #Metrum differential equation solver for pharmacometrics

# Define mrgsolve model code
  code <- '
$INIT  // Initial Conditions for Compartments
  CMT1 = 0,  // Central Compartment
  AUC = 0,  // Area under the curve

$SET  // Set Differential Equation Solver Options			
  atol      =  1e-8, rtol = 1e-8
  maxsteps  =  100000

$PARAM  // Population parameters
  POPCL = 0.179,  // Clearance (L/h/kg)
  POPV1 = 0.723,  // Central Volume (L/kg)
  
  // Default Covariate Values for Simulation (allocated in creation of population)
  BW = 70,  // 70kg
  HEIGHT = 152.4,  // 152.4cm
  SEX = 0,  // Male
  INFDUR = 1,  //Infusion Duration

  // Residual Variability
  ERR_PRO = 0.08602325,  // proportional error (fraction)
  ERR_ADD = 0.2389561,  // additive error (mg/L)

  // Default ETA Values for Simulation
  // Allocated in population so set to zero
  ETA1 = 0,  // PPVCL
  ETA2 = 0,  // PPVV1

  // Default EPS values for simulation
  // Allocated in population so set to zero
  EPS1 = 0,  // RUVPRO
  EPS2 = 0,  // RUVADD

$OMEGA   // Population Parameter Variability		
  block = TRUE
  0.0388  // PPVCL
  0.0240 0.0243  // covariance-CL,V1 + PPVV1

$SIGMA   // Residual Unexplained Variability	
  block = FALSE
  1  // Proportional error
  1  // Additive error
         
$MAIN    // Individual Parameter Values
  if((SEX == 0)&(HEIGHT > 152.4)) double IBW = 50 + 2.3*((HEIGHT - 152.4)/2.54);
  if((SEX == 0)&(HEIGHT <= 152.4)) IBW = 50;
  if((SEX == 1)&(HEIGHT > 152.4)) IBW = 45.5 + 2.3*((HEIGHT - 152.4)/2.54);
  if((SEX == 1)&(HEIGHT <= 152.4)) IBW = 45.5;
  
  if(BW > IBW) double AIBW = 0.25*(BW - IBW) + IBW;
  if(BW <= IBW) AIBW = BW;
  
  if(BW <= IBW) double CL = POPCL*BW*exp(ETA1);  // *exp(ETA(1));
  if(BW <= IBW) double V1 = POPV1*BW*exp(ETA2);  // *exp(ETA(2));
  
  if(BW > IBW) CL = POPCL*AIBW*exp(ETA1);  // *exp(ETA(1));
  if(BW > IBW) V1 = POPV1*AIBW*exp(ETA2);  // *exp(ETA(2));
  
  D_CMT1    =  INFDUR;       // Infusion duration

$ODE     // Differential Equations
  double C1 =  CMT1/V1;
  
  dxdt_CMT1 = -C1*CL;
  dxdt_AUC  =  C1;

$TABLE	 // Determines Values and Includes in Output	
  double IPRE = C1;
  double DV = IPRE*(1 + EPS1*ERR_PRO) + EPS2*ERR_ADD;  // observed concentration
  // double DV = IPRE*(1 + EPS(1)*ERR_PRO) + EPS(2)*ERR_ADD;  // debug

$CAPTURE 
  IPRE DV CL V1 ETA1 ETA2 BW HEIGHT SEX INFDUR IBW AIBW
'

# ------------------------------------------------------------------------------
# Compile the model code
  mod <- mrgsolve::mcode("BusulfanPopPKModel", code)