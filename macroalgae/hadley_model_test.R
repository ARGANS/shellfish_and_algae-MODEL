### Hadley et al., 2015 Model (DOI 10.1007/s10811-014-0370-y)

#### load libraries

library('deSolve')
library('rootSolve')
library('bvpSolve')
library('deTestSet')
library('ReacTran')
library('simecol')
library('ggplot2')

#### create harvesting data frame ######
EST <- 30 #establishment time /d
HAR <- 14 #harvesting frequency after establishment / d
HF <- 0.25 #harvest this fraction


hadley_harvest_function<-function(EST=EST,HAR=HAR,HF=HF,name,run_length) {
  num<-(run_length - EST)/HAR
  harvesttimes<-c(0,seq(from=EST,by=HAR,length.out=num))
  len<-length(harvesttimes)

  event<-data.frame(
    var=rep(name,length.out=len),
    time=harvesttimes,
    value=rep(1-HF,length.out=len),
    method=rep('multiply',length.out=len)
  )

  event


}


harvest_regime<-rbind(
  harvest_function(EST,HAR,HF,name='N_s'),
  harvest_function(EST,HAR,HF,name='N_f')
  )

#### model boundary forcings ###########################################

times <- 1:3650

### solar insolation  #umol photons m-2 s-1

PAR_mean <- 600
PAR_magn <- 400

PAR <- PAR_mean + (PAR_magn*sin(2*pi*(times)/365))
PAR <- approxfun(x = times,y = PAR, method='linear', rule=2)

### temperature  # celcius
SST_mean <- 15
SST_magn <- 3

SST <- SST_mean + (SST_magn*sin(2*pi*(times)/365))
SST <- approxfun(x = times,y = SST, method='linear', rule=2)

### NH4 # mg N m-3
NH4_mean <- 25
NH4_magn <- 10

NH4_ref <- NH4_mean + (NH4_magn*sin(2*pi*(times+180)/365))
NH4_ref <- approxfun(x = times,y = NH4_ref, method='linear', rule=2)

### NO3 # mg N m-3
NO3_mean <- 35
NO3_magn <- 10

NO3_ref <- NO3_mean + (NO3_magn*sin(2*pi*(times+180)/365))
NO3_ref <- approxfun(x = times,y = NO3_ref, method='linear', rule=2)

### Detritus # mg N m-3
D_ref  <- 0.1  #simply set to constant for now



parms <- c(
  mu      = 0.33,    # maximum growth rate           / 1/d
  V_NH4   = 60,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
  V_NO3   = 25,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
  K_NH4   = 700,     # Half saturation constant NH4  / mg(N)m-3
  K_NO3   = 300,     #   "      "          "    NO3  / mg(N)m-3
  Q_max   = 70,      # max. internal nitrogen        / mg(N) g-1 (dw)
  Q_min   = 14,      # min.    "        "            / mg(N) g-1 (dw)
  K_c     = 7,       # Half growth constant          / mg(N) g-1 (dw)
  T_O     = 12,      # optimum growth temperature    / oC
  T_r     = 1,       # range of optimum temperature  / oC
  I_s     = 277,     # saturation irradiance         / umol photons m-2 s-1
  a_cs    = 0.00036, # nitrogen-specific shading     / m2 mg-1 (N)
  d_m     = 0.003,   # mortality rate                / 1/d
  h_MA    = 0.2,     # height of seaweed             / m
  r_L     = 0.2,     # remineralisation rate         / 1/d
  r_N     = 0.1,     # nitrification rate            / 1/d
  K_d     = 0.1,     # light attenuation coefficient / 1/m
  A_farm  = 1,       # farm area                     / m2
  z       = 3,       # cultivation depth             / m
  F_in    = 3     # flow rate into farm           / m3 d-1
)

######################################################
##   THE MODEL EQUATIONS ###########
######################################################
Hadley_model <- function(t, y, parms) {

  with(as.list(c(y, parms)), {
    lambda_R    <- F_in/(A_farm*z)                                             # refresh rate of farm
    V_MA        <- h_MA * A_farm                                             # volume of macroalgae
    Q           <- Q_min*(1+(N_s/N_f))                                       # Internal nutrient quota of macroalgae
    B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
    K_MA        <- N_f*a_cs*(max(h_MA/z,1)*(min(h_MA,z)^(-1)))               # light attenuation due to macroalgae
    K           <- K_MA + K_d                                                # total light attenuation due to water and algae
    E_z         <- PAR(t)*exp(-K*z)                                          # Irradiance at top of macroalgal canopy
    g_Q         <- (Q-Q_min)/(Q-K_c)                                         # Growth limitation due to internal nutrient reserves
    g_T         <- 1/(1+exp(-(SST(t)-T_O)/T_r))                              # Growth limitation due to temperature
    g_E         <- (exp(1)/(K*h_MA))*(exp(-(E_z*exp(-K*h_MA))/I_s)-exp(-(E_z/I_s)))# Growth limitation due to light
    mu_g_EQT    <- mu*g_E*g_Q*g_T                                            # Growth function for macroalgae
    f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
    f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.
    dNH4        <- lambda_R*(NH4_ref(t)-NH4)-(f_NH4*B)+(r_L*D)-(r_N*NH4)+(d_m*N_s)     # change in NH4 with time  - eq 3 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
    dNO3        <- lambda_R*(NO3_ref(t)-NO3)-(f_NO3*B)+(r_N*NH4)             # change in NO3 with time  - eq 4 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)

    dN_s        <- (f_NH4+f_NO3)*B-mu_g_EQT*N_s-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
    dN_f        <- mu_g_EQT*N_s-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
    dD          <- lambda_R*(D_ref-D) + d_m*N_f - r_L*D                      # change in detritus with time
    list(c(dNH4, dNO3,dN_s,dN_f,dD),
         lambda_R  = lambda_R,
         V_MA      = V_MA,
         Q         = Q,
         B         = B,
         K_MA      = K_MA,
         K         = K,
         E_z       = E_z,
         g_Q       = g_Q,
         g_T       = g_T,
         g_E       = g_E,
         mu_g_EQT  = mu_g_EQT,
         f_NH4     = f_NH4,
         f_NO3     = f_NO3
         )
  })
}

y0   <- c(NH4=NH4_ref(1),NO3=NO3_ref(1),N_s=0.01,N_f=0.01,D=0.1)
start_time<-Sys.time()
 Out_Harvets <- ode(times = times, func = Hadley_model, y = y0, parms = parms, event=list(data=harvest_regime))
 Out <- ode(times = times, func = Hadley_model, y = y0, parms = parms)
 end_time<-Sys.time()

print(end_time-start_time)

plot(Out)



