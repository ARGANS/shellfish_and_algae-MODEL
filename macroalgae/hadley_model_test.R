### Hadley et al., 2015 Model (DOI 10.1007/s10811-014-0370-y)

#### load libraries

library('deSolve')
#library('rootSolve')
#library('bvpSolve')
#library('deTestSet')
#library('ReacTran')
#library('simecol')
#library('ggplot2')


boundary_forcings<-function(input_data){
  # This function takes input environmental data for an individual grid square
  # as a data frame of daily average values and generates the necessary functions
  # to return interpolated values at any time point within the range
  # (i.e. at whatever temporal resolution is required  by the ODE solver).
  # The following columns are expected in the data frame. If they are absent a
  # warning is generated and some default data substituted
  # time (in days, e.g. 1:365),
  # PAR (PAR incident at the sea surface in umol photons m-2 s-2),
  # SST (Sea surface temperature in celcius),
  # NO3 (nitrate concentration - SML average - in mg N m-3),
  # NH4 (ammonium concentration - SML average - in mg N m-3),
  # K_d (turbidity dependent light attenuation coefficient - in m-1 - derived from turbidity or SPM, satellite or modelled)
  # F_in (net horizontal flow rate into farm cross section A_xz (m/d))
  # h_z_SML (depth of SML in m)
  # t_z (vertical turnover of SML in d-1)

  maf<-make_approx_fun<-function(param){
    # given a column name construct an approxfun for that parameter with time
    assign(param,approxfun(x=input_data$time,y=input_data[,param], method='linear',rule=2),envir=.GlobalEnv)
  }

  for (name in names(input_data)) {
    maf(name)
  }

  return(invisible(0))


}




#### create harvesting data frame ######
EST <- 30 #establishment time /d
HAR <- 14 #harvesting frequency after establishment / d
HF <- 0.25 #harvest this fraction


hadley_timed_harvest_function<-function(EST=EST,HAR=HAR,HF=HF,name,run_length) {
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







#### model boundary forcings ###########################################

times <- 1:3650

make_seasonal_cycle<-function(times,mean, magnitude,offset){
  # make a sine wave seasonal cycle of a parameter given mean values,
  # magnitude (amplitude*0.5) and offset in days to time of mean value
  # as per Hadley et al 2015
  mean+(magnitude*sin(2*pi*(times+offset)/365))
}



### solar insolation  #umol phot0ons m-2 s-1

#PAR_mean <- 600
#PAR_magn <- 400

#PAR <- PAR_mean + (PAR_magn*sin(2*pi*(times)/365))
#PAR <- approxfun(x = times,y = PAR, method='linear', rule=2)

### temperature  # celcius
#SST_mean <- 15
#SST_magn <- 3

#SST <- SST_mean + (SST_magn*sin(2*pi*(times)/365))
#SST <- approxfun(x = times,y = SST, method='linear', rule=2)

### NH4 # mg N m-3
#NH4_mean <- 25
#NH4_magn <- 10
#
#NH4_ref <- NH4_mean + (NH4_magn*sin(2*pi*(times+180)/365))
#NH4_ref <- approxfun(x = times,y = NH4_ref, method='linear', rule=2)

### NO3 # mg N m-3
#NO3_mean <- 35
#NO3_magn <- 10

#NO3_ref <- NO3_mean + (NO3_magn*sin(2*pi*(times+180)/365))
#NO3_ref <- approxfun(x = times,y = NO3_ref, method='linear', rule=2)

### Detritus # mg N m-3
#D_ref  <- 0.1  #simply set to constant for now


input_hadley <- data.frame(
  time   = times,
  PAR    = make_seasonal_cycle(times,600,400,0),
  SST    = make_seasonal_cycle(times,15,3,0),
  NH4_ref = make_seasonal_cycle(times,25,20,180),
  NO3_ref = make_seasonal_cycle(times,35,10,180),
  D_ref = 0.1
)

boundary_forcings(input_hadley)


parms_ulva <- c(
  mu      = 0.45,    # maximum growth rate           / 1/d
  V_NH4   = 124,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
  V_NO3   = 39,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
  K_NH4   = 700,     # Half saturation constant NH4  / mg(N)m-3
  K_NO3   = 70,     #   "      "          "    NO3  / mg(N)m-3
  Q_max   = 42,      # max. internal nitrogen        / mg(N) g-1 (dw)
  Q_min   = 13,      # min.    "        "            / mg(N) g-1 (dw)
  K_c     = 7,       # Half growth constant          / mg(N) g-1 (dw)
  T_O     = 12,      # optimum growth temperature    / oC
  T_min     = 1,       # min temperature for growth  / oC
  T_max    = 25,     # max temperature for growth     / oC
  T_r       = 1,
  I_s     = 200,     # saturation irradiance         / umol photons m-2 s-1
  a_cs    = 0.00033, # nitrogen-specific shading     / m2 mg-1 (N)
  d_m     = 0.003,   # mortality rate                / 1/d
  h_MA    = 0.2,     # height of seaweed             / m
  r_L     = 0.2,     # remineralisation rate         / 1/d
  r_N     = 0.1     # nitrification rate            / 1/d
)

parms_porphyra <- c(
  mu      = 0.33,    # maximum growth rate           / 1/d
  V_NH4   = 60,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
  V_NO3   = 25,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
  K_NH4   = 700,     # Half saturation constant NH4  / mg(N)m-3
  K_NO3   = 300,     #   "      "          "    NO3  / mg(N)m-3
  Q_max   = 70,      # max. internal nitrogen        / mg(N) g-1 (dw)
  Q_min   = 14,      # min.    "        "            / mg(N) g-1 (dw)
  K_c     = 7,       # Half growth constant          / mg(N) g-1 (dw)
  T_O     = 12,      # optimum growth temperature    / oC
  T_min     = 1,       # min temperature for growth  / oC
  T_max    = 25,     # max temperature for growth     / oC
  T_r     = 1,
  I_s     = 277,     # saturation irradiance         / umol photons m-2 s-1
  a_cs    = 0.00036, # nitrogen-specific shading     / m2 mg-1 (N)
  d_m     = 0.003,   # mortality rate                / 1/d
  h_MA    = 0.2,     # height of seaweed             / m
  r_L     = 0.2,     # remineralisation rate         / 1/d
  r_N     = 0.1     # nitrification rate            / 1/d
)

parms_farm<-c(
  K_d     = 0.1,     # light attenuation coefficient / 1/m
  A_farm  = 1,       # farm area                     / m2
  z       = 3,       # cultivation depth             / m
  F_in    = 0.75,     # flow rate into farm           / m3 d-1
  N_farm  = 0# additional ammonium input to farm e.g. from salmon mg/N/m-3
)

print('setup done')
######################################################
##   THE MODEL EQUATIONS ###########
######################################################
Hadley_model <- function(t, y, parms,...) {

  with(as.list(c(y, parms)), {
    length_scaler<-1
    #length_scaler<-max(h_MA/z,1)
    #length_scaler<-min(h_MA/z,1)

    lambda_R    <- F_in/(A_farm*z)                                             # refresh rate of farm
    V_MA        <- h_MA * A_farm                                             # volume of macroalgae
    Q           <- Q_min*(1+(N_s/N_f))                                       # Internal nutrient quota of macroalgae
    B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
    #K_MA        <- N_f*a_cs*max(h_MA/z,1)*(min(h_MA,z)^(-1))               # light attenuation due to macroalgae
    #K_MA        <- N_f*a_cs                                                  # height scaling removed due to making no sense and not being relevant to fixed heigh MA model
    #K           <- K_MA + K_d                                                # total light attenuation due to water and algae
    #E_z         <- PAR(t)*exp(-K*z)                                          # Irradiance at top of macroalgal canopy
    g_Q         <- (Q-Q_min)/(Q-K_c)                                         # Growth limitation due to internal nutrient reserves
    #g_T         <- 1/(1+exp(-(SST(t)-T_O)/T_r))                              # Growth limitation due to temperature
    #g_E         <- (PAR(t)/(K*h_MA))*(exp(-(E_z*exp(-K*h_MA))/I_s)-exp(-(E_z/I_s)))# Growth limitation due to light

   # new light limitation - from Zollman et al 2021
    I_top<-PAR(t)*exp(-K_d*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
    I_av <-(I_top/(K_d*z+N_f*a_cs))*(1-exp(-(K_d*z+N_f*a_cs)))          # calclulate the average irradiance throughout theheight of the farm
    g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function

  # new temperature scheme introduced to give inhibition of growth at temperatureas above optimum - from martin and marques 2002

    if(SST(t)>T_O)
      {T_x<-T_max}
    else
      {T_x<-T_min}
    g_T<-exp(-2.3*((SST(t)-T_O)/(T_min-T_O))^2)

    mu_g_EQT    <- mu*g_E*g_Q*g_T                                            # Growth function for macroalgae
    f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
    f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.

    dNH4        <- lambda_R*(NH4_ref(t)-NH4)-(f_NH4*B*length_scaler)+(r_L*D)-(r_N*NH4)+(d_m*N_s)+N_farm    # change in NH4 with time  - eq 3 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
    dNO3        <- lambda_R*(NO3_ref(t)-NO3)-(f_NO3*B*length_scaler)+(r_N*NH4)             # change in NO3 with time  - eq 4 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)

    dN_s        <- (f_NH4+f_NO3)*B*length_scaler-mu_g_EQT*N_s-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley


    dN_f        <- mu_g_EQT*N_s-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
    dD          <- lambda_R*(D_ref(t)-D) + d_m*N_f - r_L*D                      # change in detritus with time
    list(c(dNH4, dNO3,dN_s,dN_f,dD),
         lambda_R  = lambda_R,
         V_MA      = V_MA,
         Q         = Q,
         B         = B,
         #K_MA      = K_MA,
         #K         = K,
         #E_z       = E_z,
         I_av     = I_av,
         g_Q       = g_Q,
         g_T       = g_T,
         g_E       = g_E,
         mu_g_EQT  = mu_g_EQT,
         f_NH4     = f_NH4,
         f_NO3     = f_NO3
         )
  })
}


Hadley_model_as_published <- function(t, y, parms) {

  with(as.list(c(y, parms)), {
    length_scaler<-max(h_MA/z,1)
    lambda_R    <- F_in/(A_farm*z)                                             # refresh rate of farm
    V_MA        <- h_MA * A_farm                                             # volume of macroalgae
    Q           <- Q_min*(1+(N_s/N_f))                                       # Internal nutrient quota of macroalgae
    B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
    K_MA        <- N_f*a_cs*length_scaler/(min(h_MA,z))               # light attenuation due to macroalgae
    K           <- K_MA + K_d                                                # total light attenuation due to water and algae
    E_z         <- PAR(t)*exp(-K*z)                                          # Irradiance at top of macroalgal canopy
    g_Q         <- (Q-Q_min)/(Q-K_c)                                         # Growth limitation due to internal nutrient reserves
    g_T         <- 1/(1+exp(-(SST(t)-T_O)/T_r))                              # Growth limitation due to temperature

    g_E         <- (exp(1)/(K*h_MA))*(exp(-(E_z*exp(-K*h_MA))/I_s)-exp(-(E_z/I_s)))# Growth limitation due to light
    mu_g_EQT    <- mu*g_E*g_Q*g_T                                            # Growth function for macroalgae
    f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
    f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.

    dNH4        <- lambda_R*(NH4_ref(t)-NH4)-(f_NH4*B*length_scaler)+(r_L*D)-(r_N*NH4)+N_farm     # change in NH4 with time  - eq 3 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
    dNO3        <- lambda_R*(NO3_ref(t)-NO3)-(f_NO3*B*length_scaler)+(r_N*NH4)             # change in NO3 with time  - eq 4 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)

    dN_s        <- (f_NH4+f_NO3)*B*length_scaler-mu_g_EQT*N_s-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley


    dN_f        <- mu_g_EQT*N_s-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
    dD          <- lambda_R*(D_ref(t)-D) + d_m*N_f - r_L*D                      # change in detritus with time
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
print('model done')
y0   <- c(NH4=NH4_ref(1),NO3=NO3_ref(1),N_s=0.1,N_f=0.1,D=0.1)
#start_time<-Sys.time()
#Out_Harvets <- ode(times = times, func = Hadley_model, y = y0, parms = parms, event=list(data=harvest_regime))
#Out <- ode(times = times, func = Hadley_model, y = y0, parms = c(parms_porphyra,parms_farm))
# end_time<-Sys.time()

#print(end_time-start_time)

#plot(Out)

#### investigating Hadley model functions


light_limitation<-function(N_f=30000,
                           a_cs=parms_porphyra['a_cs'],
                           I_s=parms_porphyra['I_s'],
                           h_MA=parms_porphyra['h_MA'],
                           K_d=parms_farm['K_d'],
                           z=parms_farm['z'],
                           PAR.=PAR(91)){


    #print(paste("N_f:",N_f))
    #print(paste("a_cs:",a_cs))
    #print(paste("z:",z))
    #print(paste("h_MA:",h_MA))
    #print(paste("I_s:", I_s))
    #print(paste("K_d:", K_d))
    #print(paste("PAR:", PAR.))


    K_MA        <- N_f*a_cs*max(h_MA/z,1)/(min(h_MA,z))

    #print(paste("K_MA:",K_MA))

    K           <- K_MA + K_d

    #print(paste("K", K))

    E_z         <- PAR.*exp(-K*z)

    #print(paste("E_z:", E_z))

    g_e<-(exp(1)/(K*h_MA))*(exp(-(E_z*exp(-K*h_MA))/I_s)-exp(-(E_z/I_s)))# Growth limitation due to light

    #print(paste("g_e:", g_e))
    g_e

}



light_limitation_zollmann<-function(N_f=30000,
                                     a_cs=parms_porphyra['a_cs'],
                                     I_s=parms_porphyra['I_s'],
                                     h_MA=parms_porphyra['h_MA'],
                                     K_d=parms_farm['K_d'],
                                     z=parms_farm['z'],
                                     PAR.=PAR(91)){

  I_top<-PAR.*exp(-K_d*(z-h_MA))
  I_av <-(I_top/(K_d*h_MA+N_f*a_cs))*(1-exp(-(K_d*h_MA+N_f*a_cs)))
  g_E <- I_av/((I_s+I_av))
  g_E

}


## =============================================================================
# harvest when light limitation is approached
## =============================================================================

rootfunc  <- function(t,y,parms,...){
  with(as.list(c(y, parms)), {
    I_top<-PAR(t)*exp(-K_d*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
    I_av <-(I_top/(K_d*z+N_f*a_cs))*(1-exp(-(K_d*z+N_f*a_cs)))          # calclulate the average irradiance throughout theheight of the farm
    g_E <- I_av/(((I_s)+I_av))

   return(g_E-harvest_threshold)  #when g_e reaches harvest threshold, return 0
  })
}

eventfunc <- function(t, y, parms,...){

  c(y[1],y[2],y[3]*(1-parms['harvest_fraction']),y[4]*(1-parms['harvest_fraction']),y[5])

  #y['N_f']y['N_f']*0.25#(1-harvest_fraction)
  #y['N_s']<-y['N_s']*0.25#(1-harvest_fraction)

}

parms_harvest<-c(
  harvest_threshold=0.2,
  harvest_fraction=0.75
)

#lim_harvest <- ode(times = times, func = Hadley_model,  y = y0, parms = c(parms_porphyra,parms_farm,parms_harvest), events=list(func=eventfunc, root=TRUE),rootfun=rootfunc)
