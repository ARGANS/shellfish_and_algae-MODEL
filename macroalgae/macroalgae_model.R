### Seaweed model for EMFF shellfish and algae project
## Based largely on Hadley et al., 2015 Model (DOI 10.1007/s10811-014-0370-y)
## Full details of scientific rationale can be found in associated ATBD
# Authors
# M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,
# Gilbert Langlois (ARGANS, Brest, France) glanglois@argans.eu


# this model contains only the code for building boundary forcing functions and
# the model equations. Model run is controlled by run_MA.R


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

  for (name in names(input_data)) {
    maf(name)
  }

  return(invisible(0))


}


######################################################
##   THE MODEL EQUATIONS ###########
######################################################
MA_model <- function(t, y, parms) {
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
    dNH4        <- lambda_R*(NH4_in(t)-NH4)-(f_NH4*B)+(r_L*D)-(r_N*NH4)+(d_m*N_s)     # change in NH4 with time  - eq 3 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
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




