############################
### Run control for macroalgal model
# Authors
# M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,...


# This script allows manual execution of the model by reading in environmental focings and parameter sets


# we need to think about how the outputs data is collected / analysed and returned to the platform
# how is the scripting set up to run the model automatically? - this script could be made executable and take config file and run name as arguments, but then we're creating a lot of temporary files - can we pass the data in straight from python?


#load the model
source('macroalgae_model.R')




# construct dummy input data
    ### solar insolation  #umol photons m-2 s-1
    PAR_mean <- 600
    PAR_magn <- 400
    ### temperature  # celcius
    SST_mean <- 15
    SST_magn <- 3
    ### NH4 # mg N m-3
       #nb 1 mgN m-3 = 1 ug N l-1
       #             = 1/14 umol l-1
       #             = 0.07 uM
       #             = 70 nM
       # so 15 to 35 mg m-3 is ~1 to 2.5 uM
    NH4_mean <- 5
    NH4_magn <- 3

    ### NO3 # mg N m-3
       # 25 to 45 mg N m-3 is 1.8 to 3.2 uM
    NO3_mean <- 60
    NO3_magn <- 50
    
    PO4_mean<-2
    PO4_magn<-1.9
    
    latitude<- 54

    run_length <- time <- 1:3650
    
    
default_input<- data.frame(
  time   = time,
  PAR    = PAR_mean + (PAR_magn*sin(2*pi*(time-91)/365)),
  SST    = SST_mean + (SST_magn*sin(2*pi*(time+180)/365)),
  NH4_in = NH4_mean + (NH4_magn*sin(2*pi*(time+91)/365)),
  NO3_in = NO3_mean + (NO3_magn*sin(2*pi*(time+91)/365)),
  PO4_in = PO4_mean + (PO4_magn*sin(2*pi*(time+91)/365)),
  K_d    = 0.1,
  F_in   = 0,
  h_z_SML= 30,
  t_z    = 0.001,
  D_in      = 0.1,
  theta = setup_solar_angle(latitude,start_day=0,ndays=length(time))
)



default_parms_ulva <- c(
  mu      = 0.45,    # maximum growth rate           / 1/d
  V_NH4   = 124,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
  V_NO3   = 39,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
  K_NH4   = 700,     # Half saturation constant NH4  / mg(N)m-3
  K_NO3   = 70,     #   "      "          "    NO3  / mg(N)m-3
  Q_max   = 42,      # max. internal nitrogen        / mg(N) g-1 (dw)
  Q_min   = 13,      # min.    "        "            / mg(N) g-1 (dw)
  N_to_P  = 12,      #N:P ratio of seaweed biomass
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

default_parms_porphyra <- c(
  mu      = 0.33,    # maximum growth rate           / 1/d
  V_NH4   = 60,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
  V_NO3   = 25,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
  K_NH4   = 700,     # Half saturation constant NH4  / mg(N)m-3
  K_NO3   = 300,     #   "      "          "    NO3  / mg(N)m-3
  Q_max   = 70,      # max. internal nitrogen        / mg(N) g-1 (dw)
  Q_min   = 14,      # min.    "        "            / mg(N) g-1 (dw)
  N_to_P  = 12,      #N:P ratio of seaweed biomass
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


default_parms_farm<-c(
  x_farm   = 1,       # width of farm in flow direction    / m
  z       = 3,       # cultivation depth             / m
  #N_farm  = 0# additional ammonium input to farm e.g. from salmon mg/N/m-3
  harvest_freq = 14, #days (only used if harvest_method==1)
  harvest_gE = 0.2, #value of light-dependent growth factor at which harvest happens (only used if harvest_method==2)
  harvest_fraction = 0.75 #fraction of total biomass to harvest (only used if narvest_method != 0)
)

default_parms_run<-c(
  refresh_rate = 1, #if value is 1, farm is fully refreshed with new water each day. Otherwise calculate from horizontal and vertical flow
  harvest_method=0 #options: 0:no varvesting, 1.fixed frequency, 2. light-driven
)





default_run <- function(parms=c(default_parms_run,default_parms_farm,default_parms_ulva),input_dat=default_input){
  boundary_forcings(input_data)
  

  y0   <- c(NH4=NH4_in(1),NO3=NO3_in(1),N_s=100,N_f=100,D=0.1)
  Out <- ode(times = input_data$time, func = MA_model, y = y0, parms = parms)
  
  
  plot(Out)
}


run_model<-function(default_parms,default_input,parms,input,y0){
  input_data<-replace(default_input,names(input),input) # use default input where no value provided
  boundary_forcings(input_data) #create boundary forcings
  parms<-replace(default_parms,names(parms),parms) #parameters use defaults where no parmater provided in run config
  Out<-ode(times = input_data$time, 
           func=MA_model,
           y=y0,
           parms=parms)
  Out
}
