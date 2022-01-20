### Seaweed model for EMFF shellfish and algae project
## Based largely on Hadley et al., 2015 Model (DOI 10.1007/s10811-014-0370-y)
## Full details of scientific rationale can be found in associated ATBD
# Authors
# M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,
# Gilbert Langlois (ARGANS, Brest, France) glanglois@argans.eu

#### load libraries

library('deSolve')
library('rootSolve')
library('bvpSolve')
library('deTestSet')
library('ReacTran')
library('simecol')
library('ggplot2')

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

  maf('PAR')
  maf('SST')
  maf('NO3_in')
  maf('NH4_in')
  maf('K_d')
  maf('F_in')
  maf('h_z_SML')
  maf('t_z')
  maf('D_in')

  return(invisible(0))


}

# construct dummy input data for testing

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

  time <- 1:730

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



default_parms <- c(
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
  #K_d     = 0.1,     # light attenuation coefficient / 1/m
  # not currently used y_farm   = 1,      # length of farm unit perpendicular to flow /m
  x_farm   = 1,       # width of farm in flow direction    / m2
  z       = 3       # cultivation depth             / m
  #F_in    = 3     # flow rate into farm           / m3 d-1
)

######################################################
##   THE MODEL EQUATIONS ###########
######################################################
Hadley_model <- function(t, y, parms) {
 #model runs per unit area - the only thing that is affected by size of the farm is the relationship
  # between flow in and width of farm and height of algae and vertical turnover i.e. nutrient resupply
  with(as.list(c(y, parms)), {
    lambda_RF    <- min(1,(F_in(t)/x_farm)) # horizontal flow refresh rate
    lambda_Rz    <-(h_z_SML(t)*t_z(t)/h_MA)*(1-lambda_RF) #vertical turnover refresh rate (in proportion of farm not refreshed by horizontal flow)
    lambda_R    <- min(1,lambda_RF+lambda_Rz) # refresh rate of farm considering both horizontal flow and vertical turnover within SML
    Q           <- Q_min*(1+(N_s/N_f))                                       # Internal nutrient quota of macroalgae
    B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
    K_MA        <- N_f*a_cs*(max(h_MA/z,1)*(min(h_MA,z)^(-1)))               # light attenuation due to macroalgae
    K           <- K_MA + K_d(t)                                                # total light attenuation due to water and algae
    E_z         <- PAR(t)*exp(-K*z)                                          # Irradiance at top of macroalgal canopy
    g_Q         <- (Q-Q_min)/(Q-K_c)                                         # Growth limitation due to internal nutrient reserves
    g_T         <- 1/(1+exp(-(SST(t)-T_O)/T_r))                              # Growth limitation due to temperature
    g_E         <- (exp(1)/(K*h_MA))*(exp(-(E_z*exp(-K*h_MA))/I_s)-exp(-(E_z/I_s)))# Growth limitation due to light
    mu_g_EQT    <- mu*g_E*g_Q*g_T                                            # Growth function for macroalgae
    f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
    f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.
    dNH4        <- lambda_R*(NH4_in(t)-NH4)-(f_NH4*B)+(r_L*D(t))-(r_N*NH4)+(d_m*N_s)     # change in NH4 with time  - eq 3 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
    dNO3        <- lambda_R*(NO3_in(t)-NO3)-(f_NO3*B)+(r_N*NH4)             # change in NO3 with time  - eq 4 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
    dN_s        <- (f_NH4+f_NO3)*B-mu_g_EQT*N_s-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
    dN_f        <- mu_g_EQT*N_s-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
    dD          <- lambda_R*(D_in(t)-D + d_m*N_f - r_L*D)                  # change in detritus with time
    list(c(dNH4, dNO3,dN_s,dN_f,dD),
         lambda_R  = lambda_R,
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



test_run <- function(parms,input_data){
  boundary_forcings(input_data)

  y0   <- c(NH4=NH4_in(1),NO3=NO3_in(1),N_s=0.01,N_f=0.01,D=0.1)
  start_time<-Sys.time()
  #Out_Harvets <- ode(times = input_data$time, func = Hadley_model, y = y0, parms = parms, event=list(data=harvest_regime))
  Out <- ode(times = input_data$time, func = Hadley_model, y = y0, parms = parms)
  end_time<-Sys.time()

  print(end_time-start_time)

  plot(Out)
}




