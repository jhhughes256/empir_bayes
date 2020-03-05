# Salinger 2010, Busulfan Population Pharmacokinetic Model

# ------------------------------------------------------------------------------
# Define the model parameters and equations

setwd("D:/Busulfan")
library(dplyr)	    #New plyr - required for mrgsolve
library(mrgsolve)	  #Metrum differential equation solver for pharmacometrics

code <- '

$INIT    // Initial Conditions for Compartments
         CMT1      =  0,         // Central Compartment
         AUC       =  0,         // Area under the curve



$SET     // Set Differential Equation Solver Options			
         atol      =  1e-8, rtol = 1e-8
         maxsteps  =  100000


$PARAM   // Population parameters
         POPCL     =  0.179,	    // Clearance (L/h/kg)
         POPV1     =  0.723,	    // Central Volume (L/kg)



         // Default Covariate Values for Simulation (allocated in creation of population)
         BW        = 70,         // 70kg
         HEIGHT    = 152.4,      // 152.4cm
         SEX       = 0,          // Male
         INFDUR    = 1,          //Infusion Duration


         // Default ETA Values for Simulation (allocated in creation of population)
         PPVCL      =  0,         // PPVCL
         PPVV1     =  0          // PPVV1


$OMEGA   // Population Parameter Variability		
         name      = "PPV"
         block     =  TRUE
         labels    =  s(ETA1, ETA2)
         0.0388                // PPVCL
         0.0240   0.0243       // covariance-CL,V1 + PPVV1


$SIGMA   // Residual Unexplained Variability	
         block     =  FALSE
         labels    =  s(ERRPROP,ERRADD)
         0.0074	               // Proportional error
         0.0571                // Additive error
         


$MAIN    // Individual Parameter Values

         if((SEX == 0)&(HEIGHT > 152.4)) double IBW = 50 + 2.3*((HEIGHT - 152.4)/2.54);
         if((SEX == 0)&(HEIGHT <= 152.4))  IBW = 50;
         if((SEX == 1)&(HEIGHT > 152.4))   IBW = 45.5 + 2.3*((HEIGHT - 152.4)/2.54);
         if((SEX == 1)&(HEIGHT <= 152.4))  IBW = 45.5;
         

         if(BW > IBW)  double AIBW = 0.25*(BW - IBW) + IBW;
         if(BW <= IBW)   AIBW = BW;

         if(BW <= IBW) double CL =  POPCL*BW*exp(PPVCL);
         if(BW <= IBW) double V1 =  POPV1*BW*exp(PPVV1);

         if(BW > IBW) CL =  POPCL*AIBW*exp(PPVCL);
         if(BW > IBW)  V1 =  POPV1*AIBW*exp(PPVV1);

         D_CMT1    =  INFDUR;       // Infusion duration


$ODE     // Differential Equations
         double C1 =  CMT1/V1;

         dxdt_CMT1 = -C1*CL;
         dxdt_AUC  =  C1;

         

$TABLE	 //Determines Values and Includes in Output	
        table(IPRE) = C1;
        table(DV)   = table(IPRE)*(1+ERRPROP)+ERRADD;

$CAPTURE 
         CL V1 ETA1 ETA2 BW HEIGHT SEX INFDUR IBW AIBW

'

# ------------------------------------------------------------------------------
# Compile the model code
mod <- mcode("BusulfanPopPKModel",code)