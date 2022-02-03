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
library('tidyverse')

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
  # theta (anlge of solar incidence degrees)

  maf<-make_approx_fun<-function(param){
    # given a column name construct an approxfun for that parameter with time
    assign(param,approxfun(x=input_data$time,y=input_data[,param], method='linear',rule=2),envir=.GlobalEnv)
  }

  for (name in names(input_data)) {
    maf(name)
  }

  return(invisible(0))


}


setup_solar_angle<-function(latitude, start_day=0, ndays){
  # Calculates max solar incidence angle theta for each day of year given latitude. 
  # Output is ready to be fed into boundary_forcings to generate approxfun
  # when start_day=0 output starts on 1 Jan
  declin<-23.45*cos(((360/365)*(1:ndays+10+start_day))*pi/180)
  90-(latitude+declin)
}


## =============================================================================
# harvesting at set frequency
## =============================================================================


# timed_harvest_function<-function(EST=EST,HAR=HAR,HF=HF,name,run_length) {
#   num<-(run_length - EST)/HAR
#   harvesttimes<-c(0,seq(from=EST,by=HAR,length.out=num))
#   len<-length(harvesttimes)
#   
#   event<-data.frame(
#     var=rep(name,length.out=len),
#     time=harvesttimes,
#     value=rep(1-HF,length.out=len),
#     method=rep('multiply',length.out=len)
#   )
#   
#   event
# }

## =============================================================================
# harvesting functions
## =============================================================================

harvest_limit_rootfunc  <- function(t,y,parms,...){
  with(as.list(c(y, parms)), {
    I_top<-PAR(t)*exp(-K_d(t)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
    I_av <-(I_top/(K_d(t)*z+N_f*a_cs))*(1-exp(-(K_d(t)*z+N_f*a_cs)))          # calclulate the average irradiance throughout theheight of the farm
    g_E <- I_av/(((I_s)+I_av))
    
    return(g_E-harvest_threshold)  #when g_e reaches harvest threshold, return 0
  })
}

harvest_timed_rootfunc <- function(t,y,parms,...){
  with(as.list(c(y, parms)), {
    return(((t-harvest_first)%%harvest_freq))
  })  
}

harvest_eventfunc <- function(t, y, parms,...){
  
  c(y[1],y[2],y[3]*(1-parms['harvest_fraction']),y[4]*(1-parms['harvest_fraction']),y[5],y[6]+y[4]*(parms['harvest_fraction']))

}

######################################################
##   THE MODEL EQUATIONS ###########
######################################################
MA_model <- function(t, y, parms) {
 #model runs per unit area - the only thing that is affected by size of the farm is the relationship
  # between flow in and width of farm and height of algae and vertical turnover i.e. nutrient resupply. 
  #The density of the farm determines the spacing between lines of width w_MA, with maximum density being 45% (spacing of w_MA between lines and 5% of area allowed for access / shipping transit through the farm etc)
  
  
  with(as.list(c(y, parms)), {
    ### Internal variables
      # farm refresh rate (i.e. nutrient resupply due to advection and vertical turnover)
      if(refresh_rate==1)
      {lambda_R=1}
      else {
        lambda_RF    <- min(1,(F_in(t)/x_farm)) # horizontal flow refresh rate
        lambda_Rz    <-(h_z_SML(t)*t_z(t)/h_MA) #vertical turnover refresh rate (in proportion of farm not refreshed by horizontal flow)
        lambda_R    <- min(1,lambda_RF+lambda_Rz) # refresh rate of farm considering both horizontal flow and vertical turnover within SML
      }
       # nutrient controls on growth, relation to biomass
      Q           <- Q_min*(1+(N_s/N_f))                                       # Internal nutrient quota of macroalgae
      B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
      g_Q         <- (Q-Q_min)/(Q-K_c)                                         # Growth limitation due to internal nutrient reserves
      # temperature-growth dynamics (from martin and marques 2002)
      if(SST(t)>T_O)
        {T_x<-T_max}
      else
        {T_x<-T_min}
      g_T<-exp(-2.3*((SST(t)-T_O)/(T_min-T_O))^2)
      
      # light limits to growth - adapted from Zollman et al 2021
      I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d(t)*z+N_f*a_cs))*(1-exp(-(K_d(t)*z+N_f*a_cs)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function

      mu_g_EQT    <- mu*g_E*g_Q*g_T                                            # Growth function for macroalgae
      f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
      f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.
      
      dNH4        <- lambda_R*(NH4_in(t)-NH4)-(f_NH4*B)+(r_L*D)-(r_N*NH4)+(d_m*N_s)     # change in NH4 with time  - eq 3 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
      dNO3        <- lambda_R*(NO3_in(t)-NO3)-(f_NO3*B)+(r_N*NH4)             # change in NO3 with time  - eq 4 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
      dN_s        <- (f_NH4+f_NO3)*B-min(mu_g_EQT*N_s,PO4_in(t)*N_to_P)-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
      dN_f        <- min(mu_g_EQT*N_s,PO4_in(t)*N_to_P)-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
      dD          <- lambda_R*(D_in(t)-D + d_m*N_f - r_L*D)                  # change in detritus with time
      dYield      <- 0
      list(c(dNH4, dNO3,dN_s,dN_f,dD,dYield),
         lambda_R  = lambda_R,
         Q         = Q,
         B         = B,
         g_Q       = g_Q,
         g_T       = g_T,
         g_E       = g_E,
         mu_g_EQT  = mu_g_EQT,
         f_NH4     = f_NH4,
         f_NO3     = f_NO3
         )
  })
}




