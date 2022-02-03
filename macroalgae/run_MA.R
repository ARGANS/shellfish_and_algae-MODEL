############################
### Run control for macroalgal model
# Authors
# M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,


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

    run_length <- time <- 1:730

dummy_input<- data.frame(
  time   = time,
  PAR    = PAR_mean + (PAR_magn*sin(2*pi*(time)/365)),
  SST    = SST_mean + (SST_magn*sin(2*pi*(time)/365)),
  NH4_in = NH4_mean + (NH4_magn*sin(2*pi*(time+180)/365)),
  NO3_in = NO3_mean + (NO3_magn*sin(2*pi*(time+180)/365)),
  K_d    = 0.1,
  F_in   = 0,
  h_z_SML= 30,
  t_z    = 0.001,
  D_in      = 0.1
)



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
  x_farm   = 1,       # width of farm in flow direction    / m
  z       = 3,       # cultivation depth             / m
  #N_farm  = 0# additional ammonium input to farm e.g. from salmon mg/N/m-3
  refresh_rate = 0 #if value is 1, farm is fully refreshed with new water each day. Otherwise calculate from horizontal and vertical flow
  
)





test_run <- function(parms,input_data){
  boundary_forcings(input_data)

  y0   <- c(NH4=NH4_in(1),NO3=NO3_in(1),N_s=0.01,N_f=0.01,D=0.1)
  Out <- ode(times = input_data$time, func = MA_model, y = y0, parms = parms)
  
  
  plot(Out)
}
