### Seaweed model for EMFF shellfish and algae project
## Based largely on Hadley et al., 2015 Model (DOI 10.1007/s10811-014-0370-y)
## Full details of scientific rationale can be found in associated ATBD
# Authors
# M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,
# Gilbert Langlois (ARGANS, Brest, France) glanglois@argans.eu


# REQUIRES R>4.1!

# this model contains only the code for building boundary forcing functions and
# the model equations. Model run is controlled by run_MA.R


#### load libraries ######

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

  output = list()
  for (param in names(input_data)) {
    output[[param]] = approxfun(x=input_data$time, y=input_data[,param], method='linear', rule=2)
  }
  
  return(output)

}




get_solar_angle<-function(latitude, doy){
  # Calculates max solar incidence angle theta for a day of year, given latitude 
  declin<-23.45*cos(((360/365)*(doy+10))*pi/180)
  90-(latitude+declin)
}



setup_solar_angle<-function(latitude, start_day=0, ndays){
  # Calculates max solar incidence angle theta for each day of year given latitude. 
  # Output is ready to be fed into boundary_forcings to generate approxfun
  # when start_day=0 output starts on 1 Jan
  declin<-23.45*cos(((360/365)*(1:ndays+10+start_day))*pi/180)
  90-(latitude+declin)
}



## =============================================================================
# deployment control (i.e. when to plant out the macroalgae)
## =============================================================================



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
  with(as.list(c(y,parms)), {
    #calculate yield accross whole farm volume i.e. scale up yield (y[6]) by volume of macroalgae
    V_MA<-y_farm*x_farm*density_MA*h_MA
    c(y[1],y[2],y[3]*(1-harvest_fraction),y[4]*(1-harvest_fraction),y[5],y[6]+y[4]*harvest_fraction*V_MA) 
  })


}




######################################################
##   run model 
######################################################

run_MA_model<-function(input,parameters,y0,output='df'){
  #function can be called from R or python, sets up boundary forcing functions from input data and executes the model function, returns a neat data frame of output values

  run_length<-max(input$time)
  times=seq(1,run_length,by=1)
  
  #create boundary forcings
  input_functions = boundary_forcings(input)
  parms = c(parameters, input_functions)

  if(parms['harvest_method']==0){
    #no harvesting
    Out<-ode(times = times,
             func=MA_model,
             y=y0,
             parms=parms)
    
  } else if(parms['harvest_method']==1){
    
    #harvest at set frequency
    
    
    Out_timed_harvest <- ode(times = times,
                             func = MA_model,
                             y = y0,
                             parms = parms,
                             event=list(func=harvest_eventfunc, root=TRUE),
                             rootfun=harvest_timed_rootfunc)
    Out<-Out_timed_harvest
    
  } else if(parms['harvest_method']==2){
    
    Out_limit_harvest <- ode(times = times,
                             func = MA_model,
                             y = y0,
                             parms = parms,
                             events=list(func=harvest_eventfunc, root=TRUE),
                             rootfun=harvest_limit_rootfunc)
    Out<-Out_limit_harvest
  }
  
  if(output=='df'){
    Out.<-as.data.frame(Out)
    cbind(input_data,Out.)
  } else {Out}
}


######################################################
##   THE MODEL EQUATIONS ###########
######################################################
MA_model <- function(t,y,parms) {
  #New model version diverges from Hadley model in its spatial structure. Hadley et al model runs notionally per unit volume of water, with flow through of water in volume units. In the old version of the model implemented here this was adapted for a flow-through and vertical water flux in m/d. A major shortcoming of the old version is the resupply of nutrients via the labmda term, which cannot be greater than 1 and therefore the total input of nutrient can only be enough per timestep to restore nutrient back to ambient levels. At a timestep of 1 day this is not sufficient given the significant drawdown of nutrients by the seaweed and the model is limited then by the lambda term (the way the Hadley model is structured values of labmda>1 when F_in > farm volume increase nutrient concentration above ambient which does not make sense). In this new version of the model, advection and vertical mixing lengths give an effective volume over which nutrient concentration change is calculated, thereby allowing for a realistic throughput of nutrient and output of meaningful delta-C.
  
  #input parms list ####
  
  # Farm dimensions: x_farm,y_farm,z_farm
  # MA dimensions: h_MA, w_MA, density_MA
  # Growth parameters: Q_min, K_c, T_O, T_max, T_min, V_NH4, V_NO3, K_NH4, K_N03, a_cs, I_s, mu, r_L, r_N, d_m,
  # MA properties: N_to_P
  
  #time varying inputs: K_d, theta,NO3_in, NH4_in, PO4_in, D_in
  
  
  with(as.list(c(y, parms)), {
    
    #dimensions ####
    #farm dimensions are defined by input parms x_farm, y_farm (both 1000m i.e. grid square dimensions by default) and z_farm (depth to bottom of macroalgae)
    #Macroalgae are grown on notional 'lines' running the full width of y_farm (perpendicular to notional flow in x direction). The 'spacing' of the lines is determined by input parm density_MA. The total volume occupied by macroalgae is calculated as y_farm*x_farm*density_MA*h_MA
    V_MA<-y_farm*x_farm*density_MA*h_MA
    # The effective volume of water that the V_MA is interacting with (concentration-dependent drawdown of nutrients) is 
    V_EFF<-y_farm*(x_farm+F_in(t))*(z+t_z(t))
    # therefore nutrient change in the outflow of the farm (to depth z+t_z) is scaled to V_MA/V_EFF of the local nutrient change (See dNH4 and dNO3 terms)
    
    lambda   <- min(1,(F_in(t)/x_farm)) #
    
    # nutrient controls on growth, relation to biomass ####
    Q           <- ifelse(N_f>0,Q_min*(1+(N_s/N_f)),0)                                       # Internal nutrient quota of macroalgae                                      # Internal nutrient quota of macroalgae
    B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
    g_Q         <- (Q-Q_min)/(Q-K_c)                                         # Growth limitation due to internal nutrient reserves
    # temperature-growth dynamics (from martin and marques 2002) ####
    if(SST(t)>T_O)
    {T_x<-T_max}
    else
    {T_x<-T_min}
    g_T<-exp(-2.3*((SST(t)-T_O)/(T_x-T_O))^2)
    
    # light limitation schemes####
    
    if(light_scheme==1){
      #light limits to growth including self-shading - adapted from Zollman et al 2021
      I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d(t)*h_MA+N_f*a_cs))*(1-exp(-(K_d(t)*h_MA+N_f*a_cs)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==2){
      # simple vertical light no shading
      I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d(t)*h_MA))*(1-exp(-(K_d(t)*h_MA)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==3){
      #solar angle accounted for no shading
      theta<-get_solar_angle(latitude,t)
      I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA)/sin(theta*pi/180))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d(t)*h_MA/sin(theta*pi/180)))*(1-exp(-(K_d(t)*h_MA/sin(theta*pi/180))))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==4){
      #solar angle accounted included with shading
      theta<-get_solar_angle(latitude,t)
      sine<-sin(theta*pi/180)
      I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA)/sine)
      I_av <-(I_top / ((K_d(t)*h_MA/sine)   +   (N_f*a_cs/(sine))))*(1-exp(-(  (K_d(t)*h_MA/sine)  +   (N_f*a_cs/(sine))   )))
      g_E <- I_av/(((I_s)+I_av)) 
    } 
    
    
    # growth function (uptake of nutrient from internal reserves) ####
    
    mu_g_EQT    <- mu*g_E*g_Q*g_T                                            # Growth function for macroalgae
    
    
    # uptake of nutrients from the water ####
    
    f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
    f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.
    
    # amount of nutrient removed per m3 
    NO3_removed <- (f_NO3*B)-(r_N*NH4)
    NH4_removed <- (f_NH4*B)-(r_L*D)+(r_N*NH4)-(d_m*N_s)
    
    #How much phosphate is available for growth
    PO4_tot<-PO4_in(t)*V_EFF/V_MA
    #convert N:P from mol/mol to g/g
    N_to_P<-N_to_P*14/31
    
    #model state variables ####
    dNH4        <- min(1,lambda)*(NH4_in(t)-NH4) -((f_NH4*B)-(r_L*D)+(r_N*NH4)-(d_m*N_s))*V_MA/V_EFF  # change in NH4 with time
    
    dNO3        <- min(1,lambda)*(NO3_in(t)-NO3) -((f_NO3*B)-(r_N*NH4))*V_MA/V_EFF           # change in NO3 with time 
    
    dN_s        <- (f_NH4+f_NO3)*B-min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
    
    dN_f        <- min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
    
    dD          <- min(1,lambda)*(D_in(t)-D) + (d_m*N_f - r_L*D)*V_MA/V_EFF                 # change in detritus with time
    
    
    dYield      <- 0
    
    output<-list(c(dNH4, dNO3,dN_s,dN_f,dD,dYield),
         Farm_NO3_demand = NO3_removed*V_MA,
         Farm_NH4_demand = NH4_removed*V_MA,
         V_EFF = V_EFF,
         volscale=V_EFF/V_MA,
         Q         = Q,
         B_vol  = B,
         B_line = B*h_MA*w_MA,
         g_Q       = g_Q,
         g_T       = g_T,
         g_E       = g_E,
         mu_g_EQT  = mu_g_EQT,
         f_NH4     = f_NH4,
         f_NO3     = f_NO3
    )
  }) 
}

# MA_model_old <- function(t, y, parms) {
#  #model runs per unit area - the only thing that is affected by size of the farm is the relationship
#   # between flow in and width of farm and height of algae and vertical turnover i.e. nutrient resupply. 
#   #The density of the farm determines the spacing between lines of width w_MA, with maximum density being 45% (spacing of w_MA between lines and 5% of area allowed for access / shipping transit through the farm etc)
#   
#   
#   with(as.list(c(y, parms)), {
#     ### Internal variables
#       # farm refresh rate (i.e. nutrient resupply due to advection and vertical turnover)
#       if(refresh_rate==1)
#       {lambda_R=1}
#       else {
#         lambda_RF    <- min(1,(F_in(t)/x_farm)) # horizontal flow refresh rate
#         lambda_Rz    <-(t_z(t)/h_MA) #vertical turnover refresh rate (in proportion of farm not refreshed by horizontal flow)
#         lambda_R    <- min(1,lambda_RF+(1-lambda_RF)*lambda_Rz) # refresh rate of farm considering both horizontal flow and vertical turnover within SML
#       }
#        # nutrient controls on growth, relation to biomass
#       Q           <- Q_min*(1+(N_s/N_f))                                       # Internal nutrient quota of macroalgae
#       B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
#       g_Q         <- (Q-Q_min)/(Q-K_c)                                         # Growth limitation due to internal nutrient reserves
#       # temperature-growth dynamics (from martin and marques 2002)
#       if(SST(t)>T_O)
#         {T_x<-T_max}
#       else
#         {T_x<-T_min}
#       g_T<-exp(-2.3*((SST(t)-T_O)/(T_x-T_O))^2)
#       
#       if(light_scheme==1){
#       #light limits to growth including self-shading - adapted from Zollman et al 2021
#         I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
#         I_av <-(I_top/(K_d(t)*h_MA+N_f*a_cs))*(1-exp(-(K_d(t)*h_MA+N_f*a_cs)))          # calclulate the average irradiance throughout the height of the farm
#         g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
#       } else if (light_scheme==2){
#       # simple vertical light no shading
#           I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
#           I_av <-(I_top/(K_d(t)*h_MA))*(1-exp(-(K_d(t)*h_MA)))          # calclulate the average irradiance throughout the height of the farm
#           g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
#       } else if (light_scheme==3){
#         #solar angle accounted for no shading
#           I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA)/sin(theta(t)*pi/180))                                    # calculate incident irradiance at the top of the farm
#           I_av <-(I_top/(K_d(t)*h_MA/sin(theta(t)*pi/180)))*(1-exp(-(K_d(t)*h_MA/sin(theta(t)*pi/180))))          # calclulate the average irradiance throughout the height of the farm
#           g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
#       } 
#       
#       
#       mu_g_EQT    <- mu*g_E*g_Q*g_T                                            # Growth function for macroalgae
#       f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
#       f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.
#       
#       NH4_added      <-lambda_R*(NH4_in(t)-NH4)
#       NO3_added      <-lambda_R*(NO3_in(t)-NO3)
#       
#       NO3_removed <- (f_NO3*B)-(r_N*NH4)
#       NH4_removed <- (f_NH4*B)-(r_L*D)+(r_N*NH4)-(d_m*N_s)
#       
# 
#       dNH4        <- lambda_R*(NH4_in(t)-NH4)-((f_NH4*B)+(r_L*D)-(r_N*NH4)+(d_m*N_s))/30   # change in NH4 with time  - eq 3 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
#       dNO3        <- lambda_R*(NO3_in(t)-NO3)-((f_NO3*B)+(r_N*NH4))/30            # change in NO3 with time  - eq 4 in Hadley et al (note max(h_MA/z,1) term is omitted because we assume the surface box is well mixed - need to think further about this - should we be looking across entire mixed layer - or include this in the lambda term?)
#       dN_s        <- (f_NH4+f_NO3)*B-min(mu_g_EQT*N_f,PO4_in(t)*N_to_P)-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
#       dN_f        <- min(mu_g_EQT*N_f,PO4_in(t)*N_to_P)-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
# 
# 
#       dN_s        <- (f_NH4+f_NO3)*B-min(mu_g_EQT*N_s,PO4_in(t)*N_to_P)-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
#       dN_f        <- min(mu_g_EQT*N_s,PO4_in(t)*N_to_P)-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
# 
#       dD          <- lambda_R*(D_in(t)-D + d_m*N_f - r_L*D)                  # change in detritus with time
#       dYield      <- 0
#       list(c(dNH4, dNO3,dN_s,dN_f,dD,dYield),
#          NH4_added = NH4_added,
#          NO3_added = NO3_added,
#          NO3_removed = NO3_removed,
#          NH4_removed = NH4_removed,
#          lambda_R  = lambda_R,
#          Q         = Q,
#          B         = B,
#          g_Q       = g_Q,
#          g_T       = g_T,
#          g_E       = g_E,
#          mu_g_EQT  = mu_g_EQT,
#          f_NH4     = f_NH4,
#          f_NO3     = f_NO3
#          )
#   })
# }
# 



