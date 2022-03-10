###### Run control for macroalgal model ###################
# Authors
# M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,...


# This script allows manual execution of the model by reading in environmental focings and parameter sets


#required libraries ####
library(rjson)


#load the model ####
source('macroalgae_model.R')


# dummy data (in lieu of environmental forcings) ######## -------------------------

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
    
    PO4_mean<-50
    PO4_magn<-3
    
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
  F_in   = 100,
  t_z    = 10,
  D_in   = 0
)


## parms required by macroalge_model:
#input parms list ####

# Farm dimensions: latitdue, x_farm,y_farm,z_farm
# MA dimensions: h_MA, w_MA, density_MA
# Growth parameters: Q_min, K_c, T_O, T_max, T_min, V_NH4, V_NO3, K_NH4, K_N03, a_cs, I_s, mu, r_L, r_N, d_m,
# MA properties: N_to_P

#load in all parameters from json file
allparams<-fromJSON(file='model_parameters.json')

parms_ulva<-unlist(allparams$species$ulva$parameters)
parms_saccharina<-unlist(allparams$species$saccharina$parameters)
parms_alaria<-unlist(allparams$species$alaria$parameters)

default_parms_farm<-unlist(allparams$farm$default$parameters)
default_parms_run<-unlist(allparams$run$default$options)

harvest_CCA_run<-unlist(allparams$harvest$CCA$parameters)
harvest_winter_growth_run<-unlist(allparams$harvest$Winter_growth$parameters)




# deafult params the old way ------------------------------------
# Default algae parameters: Ulva 

# default_parms_ulva <- c(
#   mu      = 0.45,    # maximum growth rate           / 1/d
#   V_NH4   = 124,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
#   V_NO3   = 39,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
#   K_NH4   = 700,     # Half saturation constant NH4  / mg(N)m-3
#   K_NO3   = 70,     #   "      "          "    NO3  / mg(N)m-3
#   Q_max   = 42,      # max. internal nitrogen        / mg(N) g-1 (dw)
#   Q_min   = 13,      # min.    "        "            / mg(N) g-1 (dw)
#   N_to_P  = 12,      #N:P ratio of seaweed biomass
#   K_c     = 7,       # Half growth constant          / mg(N) g-1 (dw)
#   T_O     = 12,      # optimum growth temperature    / oC
#   T_min     = 1,       # min temperature for growth  / oC
#   T_max    = 25,     # max temperature for growth     / oC
#   I_s     = 200,     # saturation irradiance         / umol photons m-2 s-1
#   a_cs    = 0.00033, # nitrogen-specific shading     / m2 mg-1 (N)
#   d_m     = 0.0003,   # mortality rate                / 1/d
#   h_MA    = 0.2,     # height of seaweed             / m
#   w_MA    = 1,     # width of seaweed e.g. on rope /m
#   r_L     = 0.1,     # remineralisation rate         / 1/d
#   r_N     = 0.1     # nitrification rate            / 1/d
# )


# # Default algae parameters: Alaria
# default_parms_saccharina <- c(
#   mu      = 0.06,    # maximum growth rate           / 1/d
#   V_NH4   = 100,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
#   V_NO3   = 200,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
#   K_NH4   = 11,     # Half saturation constant NH4  / mg(N)m-3
#   K_NO3   = 200,     #   "      "          "    NO3  / mg(N)m-3
#   Q_max   = 22,      # max. internal nitrogen        / mg(N) g-1 (dw)
#   Q_min   = 10,      # min.    "        "            / mg(N) g-1 (dw)
#   N_to_P  = 12,      #N:P ratio of seaweed biomass
#   K_c     = 8,       # Half growth constant          / mg(N) g-1 (dw)
#   T_O     = 12.5,      # optimum growth temperature    / oC
#   T_min     = 0,       # min temperature for growth  / oC
#   T_max    = 20,     # max temperature for growth     / oC
#    I_s     =90,     # saturation irradiance         / umol photons m-2 s-1
#   a_cs    = 0.00036, # nitrogen-specific shading     / m2 mg-1 (N)
#   d_m     = 0.0003,   # mortality rate                / 1/d
#   h_MA    = 2,     # height of seaweed             / m
#   w_MA    = 0.3,     # width of seaweed e.g. on rope /m
#   r_L     = 0.10,     # remineralisation rate         / 1/d
#   r_N     = 0.1     # nitrification rate            / 1/d
# )
# 
# 
# 
# # Default algae parameters: Porphyra 
# 
# default_parms_porphyra <- c(
#   mu      = 0.33,    # maximum growth rate           / 1/d
#   V_NH4   = 60,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
#   V_NO3   = 25,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
#   K_NH4   = 700,     # Half saturation constant NH4  / mg(N)m-3
#   K_NO3   = 100,     #   "      "          "    NO3  / mg(N)m-3
#   Q_max   = 70,      # max. internal nitrogen        / mg(N) g-1 (dw)
#   Q_min   = 14,      # min.    "        "            / mg(N) g-1 (dw)
#   N_to_P  = 12,      #N:P ratio of seaweed biomass
#   K_c     = 7,       # Half growth constant          / mg(N) g-1 (dw)
#   T_O     = 12,      # optimum growth temperature    / oC
#   T_min     = 1,       # min temperature for growth  / oC
#   T_max    = 25,     # max temperature for growth     / oC
#   I_s     = 277,     # saturation irradiance         / umol photons m-2 s-1
#   a_cs    = 0.00036, # nitrogen-specific shading     / m2 mg-1 (N)
#   d_m     = 0.003,   # mortality rate                / 1/d
#   h_MA    = 0.4,     # height of seaweed             / m
#   w_MA    = 0.2,     # width of seaweed e.g. on rope /m
#   r_L     = 0.2,     # remineralisation rate         / 1/d
#   r_N     = 0.1     # nitrification rate            / 1/d
# )



# Default Farm Parms 
# 
# default_parms_farm<-c(
#   latitude=52,
#   y_farm = 1000,       # width of farm perpendicular to flow direction    
#   density_MA = 0.4,      # fraction of farm area occupied by algae
#   x_farm = 1000,            #farm length in flow direction  
#   z       = 2,       # cultivation depth             / m
# 
#   harvest_first = 60, #days from start of run to first harvest
#   harvest_freq = 30, #days (only used if harvest_method==1)
#   harvest_threshold = 0.2, #value of light-dependent growth factor (g_E) at which harvest happens (only used if harvest_method==2)
#   harvest_fraction = 0.75 #fraction of total biomass to harvest (only used if harvest_method != 0)
# )
# 







# Default run parms 
# 
# default_parms_run<-c(
#   #refresh_rate = 1, #if value is 1, farm is fully refreshed with new water each day. Otherwise calculate from horizontal and vertical flow
#   harvest_method=0, #options: 0:no harvesting, 1.fixed frequency, 2. light-driven
#   light_scheme=4 #options 1: Zollman self-shading scheme, 2: simple vertical light no self shading, 3: solar angle light no self shading,4: selfshading with solar angle accounted for. 
#   )
# 


# test parms ulva ============
test_parms_ulva <- c(
  mu      = 0.45,    # maximum growth rate           / 1/d
  V_NH4   = 128,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
  V_NO3   = 40,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
  K_NH4   = 100,     # Half saturation constant NH4  / mg(N)m-3
  K_NO3   = 70,     #   "      "          "    NO3  / mg(N)m-3
  Q_max   = 42,      # max. internal nitrogen        / mg(N) g-1 (dw)
  Q_min   = 13,      # min.    "        "            / mg(N) g-1 (dw)
  N_to_P  = 12,      #N:P ratio of seaweed biomass
  K_c     = 7,       # Half growth constant          / mg(N) g-1 (dw)
  T_O     = 12,      # optimum growth temperature    / oC
  T_min     = 1,       # min temperature for growth  / oC
  T_max    = 25,     # max temperature for growth     / oC
  I_s     = 200,     # saturation irradiance         / umol photons m-2 s-1
  a_cs    = 0.00033, # nitrogen-specific shading     / m2 mg-1 (N)
  d_m     = 0.003,   # mortality rate                / 1/d
  h_MA    = 0.2,     # height of seaweed             / m
  w_MA    = 1,     # width of seaweed e.g. on rope /m
  r_L     = 0.1,     # remineralisation rate         / 1/d
  r_N     = 0.1     # nitrification rate            / 1/d
)

## make full defulat parameter set --------
default_parms<-c(default_parms_run,default_parms_farm,parms_ulva,harvest_winter_growth_run)

setup_run_parms<-function(default_parms.=default_parms, parms){
  #when calling rum_MA_model from R, use this function to overwrite default parameters as needed
  run_parms<-replace(default_parms.,names(parms),parms)
  run_parms
}

setup_run_input<-function(default_input.=default_input,input){
  df<-default_input.[1:nrow(input),]
  df[,colnames(input)]<-input
  df
}


# Idealised reference runs
#-------------------------------
  
# These are run using the reference_run() function below by applying different input_data as defined here for the 4 reference runs. 
  
##1 Low nutrient, low flow
  
lnlf<-c(
  NH4_mean = 2,
  NH4_magn = 1.8,
  NO3_mean = 10,
  NO3_magn = 8,
  PO4_mean=50,
  PO4_magn=30,
  F_in=10,
  t_z=5
)
  

lnhf<-c(
  NH4_mean = 2,
  NH4_magn = 1.8,
  NO3_mean = 10,
  NO3_magn = 8,
  PO4_mean=50,
  PO4_magn=30,
  F_in=1000,
  t_z=10
)
  
  
hnhf<-c(
  NH4_mean = 60,
  NH4_magn = 30,
  NO3_mean = 150,
  NO3_magn = 90,
  PO4_mean=50,
  PO4_magn=30,
  F_in=1000,
  t_z=10
)


hnlf<-c(
  NH4_mean = 60,
  NH4_magn = 30,
  NO3_mean = 150,
  NO3_magn = 90,
  PO4_mean=50,
  PO4_magn=30,
  F_in=0,
  t_z=0
)


### override default parms for reference harvest runs
parms_ref_harvest_run_winter_growth<-c(
  harvest_method=1,
  harvest_winter_growth_run
)

### override default parms for reference harvest runs
parms_ref_harvest_run_CCA<-c(
  harvest_method=2,
  harvest_CCA_run
)

# reference run function ----------------------------------------


reference_run <- function(input_data,nondefault_parms,harvest=FALSE){
  with(as.list(input_data), {
  time<-1:730
  input_frame<- data.frame(
    time   = time,
    PAR    = 400 + (300*sin(2*pi*(time-91)/365)),
    SST    = 13 + (4*sin(2*pi*(time+180)/365)),
    NH4_in = NH4_mean + (NH4_magn*sin(2*pi*(time+91)/365)),
    NO3_in = NO3_mean + (NO3_magn*sin(2*pi*(time+91)/365)),
    PO4_in = PO4_mean + (PO4_magn*sin(2*pi*(time+91)/365)),
    K_d    = 0.1,
    F_in   = F_in,
    D_in   = 0.1,
    t_z    = t_z
  )
  
  if(harvest==TRUE){
    y0   <- c(NH4=input_frame$NH4_in[1],NO3=input_frame$NO3_in[1],N_s=0,N_f=0,D=0,Yield_farm=0,Yield_per_m=0)
  }else{
    y0   <- c(NH4=input_frame$NH4_in[1],NO3=input_frame$NO3_in[1],N_s=0,N_f=100,D=0,Yield_farm=0,Yield_per_m=0)
  }
  
  parms.=setup_run_parms(parms=nondefault_parms)
  run_MA_model(input=input_frame,parameters = parms.,y0=y0,output='odeout')
  })
  
}




#lnlf_out<-reference_run(lnlf,test_parms_ulva)
#hnlf_out<-reference_run(hnlf,test_parms_ulva)
#lnhf_out<-reference_run(lnhf,test_parms_ulva)
#hnhf_out<-reference_run(hnhf,test_parms_ulva)


# lnlf_harvest_out<-reference_run(lnlf,c(test_parms_ulva,parms_ref_harvest_run),harvest=TRUE)
# hnlf_harvest_out<-reference_run(hnlf,c(test_parms_ulva,parms_ref_harvest_run),harvest=TRUE)
# lnhf_harvest_out<-reference_run(lnhf,c(test_parms_ulva,parms_ref_harvest_run),harvest=TRUE)
# hnhf_harvest_out<-reference_run(hnhf,c(test_parms_ulva,parms_ref_harvest_run_winter_growth),harvest=TRUE)







# # Main function: run_model-------------------------------------------
# 
# run_model<-function(default_parms,default_input,parms=NULL,input=NULL,y0){
#   
#   run_length<-max(input_data$time)
#   times=seq(1,run_length,by=1)
#   # use default input where no values provided
#   input_functions = boundary_forcings(input_data) #create boundary forcings
#   parms = c(parms, input_functions)
#   parms<-replace(default_parms,names(parms),parms) #parameters use defaults where no parmater provided in run config
#   if(parms['harvest_method']==0){
#     #no harvesting
#     Out<-ode(times = times,
#              func=MA_model,
#              y=y0,
#              parms=parms)
# 
#   } else if(parms['harvest_method']==1){
# 
#     #harvest at set frequency
# 
# 
#     Out_timed_harvest <- ode(times = times,
#                              func = MA_model,
#                              y = y0,
#                              parms = parms,
#                              event=list(func=harvest_eventfunc, root=TRUE),
#                              rootfun=harvest_timed_rootfunc)
#     Out<-Out_timed_harvest
# 
#   } else if(parms['harvest_method']==2){
# 
#     Out_limit_harvest <- ode(times = times,
#                              func = MA_model,
#                              y = y0,
#                              parms = parms,
#                              events=list(func=harvest_eventfunc, root=TRUE),
#                              rootfun=harvest_limit_rootfunc)
#     Out<-Out_limit_harvest
#   }
#   #cbind(input_data,Out)
#   Out
# }
