# Colin et al, 2019: Vancomycin Pharmacokinetics Throughout Life

# ------------------------------------------------------------------------------
# Define the model parameters and equations

library(plyr)
library(dplyr)
library(mrgsolve)
library(ggplot2)

code  <- '

$INIT    // Initial Conditions for Compartments
         CMT1      =  0,                         // Central Compartment
         CMT2      =  0,                         // Peripheral Compartment
         AUC       =  0,                         // Area under the curve


$SET     // Set Differential Equation Solver Options			
         atol      =  1e-8, rtol = 1e-8
         maxsteps  =  100000


$PARAM   // Population parameters
         POPCL     =  exp(1.67),                 // Clearance = 5.31 L/hr, Theta(1)
         POPV1     =  exp(3.75),                 // Central Volume = 42.9 L, Theta(2)
         POPV2     =  exp(3.73),                 // Peripheral Volume = 41.7 L, Thate(3)
         POPQ      =  exp(1.17),                 // Intercompartmental Clearance  = 3.22 L/hr, Theta(4)


         // Covariate Effects
         AGE50_CL  =  exp(4.12),                 // Age at 50% decline in CL = 61.6 years
         PMA50_CL  =  exp(-0.113),               // Age at 50% maturation in CL = 0.893 years (46.4 weeks)
         STEEP1    =  exp(1.06),                 // Hill coefficient for maturation in CL = 2.89
         STEEP2    = -exp(0.806),                // Hill coefficient for decline in CL = -2.24
         SCR_CL    =  exp(-0.436),               // Effect of SCr on CL = 0.649
         HEEL_V1   =  0.313,                     // Effect of finger-prick sampling on V1, STDY13_V1
         HEEL_Q    =  0.595,                     // Effect of finger-prick sampling on Q, STDY13_Q
         CANCER_CL =  0.293,                     // Effect of cancer on CL, STDY10_CL


         // Default Covariate Values for Simulation
         AGE       =  35,                        // Postnatal age = 35 years
         WT        =  70,                        // Body weight = 70 kg
         SCR       =  0.84,                      // Serum creatinine = 0.84 mg/dL
         HEEL      =  0,                         // Finger-prick sampling, No = 0
         CANCER    =  0,                         // Cancer, No = 0
         INFDUR    =  1,                         // Infusion Duration = 1 hour


         // Default ETA Values for Simulation
         PPVCL     = 0,                          // Population parameter variability in CL, ETA(1)                     
         PPVV1     = 0,                          // Population parameter variability in V1, ETA(2)                         
         PPVV2     = 0,                          // Population parameter variability in V2, ETA(3)                        


$OMEGA   // Population Parameter Variability		
         name      = "PPV"
         block     =  FALSE
         labels    =  s(ETA1, ETA2, ETA3)
         0.0749                                  // PPVCL = 27.4% CV, values of the NONMEM code 
         0.0721                                  // PPVV1 = 26.9% CV
         0.6720                                  // PPVV2 = 82.0% CV


$SIGMA   // Residual Unexplained Variability	
         block     =  FALSE
         labels    =  s(ERRPROP, ERRADD)
         0.046                                   // Proportional error
         1.230                                   // Additional error  


$MAIN    // Individual Parameter Values
         double Fsize   =  WT/70;                                                                  // Factor for relative body weight
         
         double PMA     =  ((AGE*365)+(40*7))/365;                                                     // Postmenstral age in years, assume gestational age = 40 weeks
         double Fmat    =  pow(PMA,STEEP1)/(pow(PMA,STEEP1)+pow(PMA50_CL,STEEP1));                       // Factor for PMA-related maturation in CL
         double Fdec    =  pow(PMA,STEEP2)/(pow(PMA,STEEP2)+pow(AGE50_CL,-STEEP2));                       // Factor for PMA-related decline in CL

         double SCRstd  =  exp(-1.22839 + log10(PMA)*0.67190 + 6.27017*exp(-3.10940*PMA));         // Reference SCr for postmenstral age
         double FSCr    =  exp(-SCR_CL*(SCR-SCRstd));                                              // Factor for SCr-related change in CL

         double CL      =  POPCL*pow((V1/POPV1),0.75)*Fmat*Fdec*FSCr*(1 + CANCER*CANCER_CL)*exp(PPVCL);
         double V1      =  POPV1*Fsize*(1-HEEL*HEEL_V1)*exp(PPVV1);
         double V2      =  POPV2*Fsize*exp(PPVV2);
         double Q       =  POPQ*pow((V2/POPV2),0.75)*(1-HEEL*HEEL_Q);

         D_CMT1         =  INFDUR;


$ODE     // Differential Equations
         double C1      =  CMT1/V1;
         double C2      =  CMT2/V2;

         dxdt_CMT1      = -C1*Q + C2*Q - C1*CL;
         dxdt_CMT2      =  C1*Q - C2*Q;
         dxdt_AUC       =  C1;     


$TABLE	 //Determines Values and Includes in Output	
         table(IPRE)    =  C1;
         table(DV)      =  table(IPRE)*((1+ERRPROP)+ERRADD);

$CAPTURE  
         CL V1 V2 Q PPVCL PPVV1 PPVV2 AGE WT SCR HEEL CANCER INFDUR PMA Fmat Fdec SCRstd FSCr
         '

# ------------------------------------------------------------------------------
# Compile the model code
mod <- mcode("VancomycinPopPKModel",code)