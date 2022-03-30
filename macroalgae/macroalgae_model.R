### Seaweed model for EMFF shellfish and algae project
## Based  on Hadley et al., 2015 Model (DOI 10.1007/s10811-014-0370-y)
## Full details of scientific rationale can be found in associated ATBD
# Authors
# M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,
# Gilbert Langlois (ARGANS, Brest, France) glanglois@argans.eu
# Quentin, Simona, Dee...


# REQUIRES R>4.1!

# to control this model using R, refer to run_MA.


#### load libraries ######

library('deSolve')
library('rootSolve')
#library('bvpSolve')
#library('deTestSet')
#library('ReacTran')
#library('simecol')
library('tidyverse')

boundary_forcings<-function(input_data, method='linear'){
  # This function takes input environmental data for an individual grid square
  # as a data frame of daily average values and generates the necessary functions
  # to return interpolated values at any time point within the range
  # (i.e. at whatever temporal resolution is required  by the ODE solver).
  # The following columns are expected in the data frame. 
  
  # time (in days, e.g. 1:365),
  # PAR (PAR incident at the sea surface in umol photons m-2 s-2),
  # SST (Sea surface temperature in celcius),
  # NO3_ext (nitrate concentration - SML average - in mg N m-3),
  # NH4_ext (ammonium concentration - SML average - in mg N m-3),
  # K_d (turbidity dependent light attenuation coefficient - in m-1 - derived from turbidity or SPM, satellite or modelled)
  # F_in (net horizontal flow rate into farm cross section A_xz (m/d))
  # t_z (vertical turnover of SML in md-1)


  output = list()
  for (param in names(input_data)) {
    output[[param]] = approxfun(x=input_data$time, y=input_data[,param], method=method, rule=2)
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
# harvesting conrol functions - when to plant out and harvest the seaweed
## =============================================================================


harvest_timed_rootfunc <- function(t,y,parms,...){
  with(as.list(c(y, parms)), {
    yroot<-c(t-deployment_day,(t-harvest_first)%%harvest_freq)
    return(yroot)
  })  
}
 
harvest_eventfunc <- function(t, y, parms,...){
  with(as.list(c(y,parms)), {
    if(t-deployment_day==0){
      #deploy some seaweed
      c(y[1],y[2],y[3],y[4]+deployment_Nf,y[5],y[6],y[7]) 
    }  else {
      #harvest; calculate yield either accross whole farm volume i.e. scale up yield (y[6]) by volume of macroalgae or per m of line (y[7]). In both cases divide by Q_min to get DW Biomass in g
      V_MA<-y_farm*x_farm*density_MA*h_MA
      c(y[1],y[2],y[3]*(1-harvest_fraction),y[4]*(1-harvest_fraction),y[5],y[6]+y[4]*harvest_fraction*V_MA/Q_min,y[7]+y[4]*harvest_fraction*h_MA*w_MA/Q_min) 
    }
  })


}




######################################################
##   run model 
######################################################

run_MA_model<-function(input,parameters,y0,output='df'){
  #function can be called from R or python, sets up boundary forcing functions from input data and executes the model function, returns a neat data frame of output values

  #create boundary forcings
  input_functions = boundary_forcings(input, "constant")
  parms = c(parameters, input_functions)

  if (!is.numeric(y0)) {
    y0_names = names(y0)
    y0 = as.numeric(y0)
    names(y0) = y0_names
  }

  if(parms['harvest_method']==0){
    #no harvesting
    Out<-ode(times = input$time,
             func=MA_model,
             y=y0,
             parms=parms)

    
  } else if(parms['harvest_method']==1){
    
    #harvest at set frequency
    Out_timed_harvest <- ode(times = input$time,
                             func = MA_model,
                             y = y0,
                             parms = parms,
                             event=list(func=harvest_eventfunc, root=TRUE),
                             rootfun=harvest_timed_rootfunc)
    Out<-Out_timed_harvest
    
  } else if(parms['harvest_method']==2){
    
    Out_limit_harvest <- ode(times = input$time,
                             func = MA_model,
                             y = y0,
                             parms = parms,
                             events=list(func=harvest_eventfunc, root=TRUE),
                             rootfun=harvest_timed_rootfunc)
    Out<-Out_limit_harvest
  }
  
  if(output=='df'){
    return(as.data.frame(Out))
  } else {Out}
}


######################################################
##   THE MODEL EQUATIONS ###########
######################################################
MA_model <- function(t,y,parms,...) {
  #New model version diverges from Hadley model in its spatial structure. Hadley et al model runs notionally per unit volume of water, with flow through of water in volume units. In the old version of the model implemented here this was adapted for a flow-through and vertical water flux in m/d. A major shortcoming of the old version is the resupply of nutrients via the labmda term, which cannot be greater than 1 and therefore the total input of nutrient can only be enough per timestep to restore nutrient back to ambient levels. At a timestep of 1 day this is not sufficient given the significant drawdown of nutrients by the seaweed and the model is limited then by the lambda term (the way the Hadley model is structured values of labmda>1 when F_in > farm volume increase nutrient concentration above ambient which does not make sense). In this new version of the model, advection and vertical mixing lengths give an effective volume over which nutrient concentration change is calculated, thereby allowing for a realistic throughput of nutrient and output of meaningful delta-C.
  
  #parameters expected: ####
  # Farm location:
  #  latitude (degrees)
  # Farm dimensions: 
  #  x_farm (m), 
  #  y_farm (m),
  #  z (m)
  # MA dimensions: 
  #  h_MA (m), 
  #  w_MA (m), 
  #  density_MA (fraction - 0 to 1)
  # Growth parameters: 
  #  Q_min (mg N per g dry weight), 
  #  K_c (mg N per g dry weight), 
  #  T_O (degrees celcius), 
  #  T_max (degrees celcius), 
  #  T_min (degrees celcius), 
  #  V_NH4 (mg N per g dry weight per day), 
  #  V_NO3 (mg N per g dry weight per day), 
  #  K_NH4 (mg N per m^3), 
  #  K_N03 (mg N per m^3), 
  #  a_cs  (m^2 per mg N), 
  #  I_s,  (micromol photons per m^2 per s)
  #  mu,  1/d
  #  r_L, 1/d
  #  r_N, 1/d
  #  d_m, 1/d
  # MA properties: 
  #  N_to_P, molar ratio (unitless),
  
  #time varying inputs: 
  # K_d, 1/m 
  # SST, degrees celcius
  # NO3_ext, mg N per m^3
  # NH4_ext, mg N per m^3
  # PO4_ext, mmol P per m^3 (i.e. uM)
  # D_ext,  mg N per m^3
  # t_z, m
  #F_in, m/d
  
  
  with(as.list(c(y, parms)), {
    
    #dimensions ####
    #farm dimensions are defined by input parms x_farm, y_farm (both 1000m i.e. grid square dimensions by default) and z_farm (depth to bottom of macroalgae)
    #Macroalgae are grown on notional 'lines' running the full width of y_farm (perpendicular to notional flow in x direction). The 'spacing' of the lines is determined by input parm density_MA. The total volume occupied by macroalgae is calculated as y_farm*x_farm*density_MA*h_MA
    V_MA<-y_farm*x_farm*density_MA*h_MA
    # The volume of water that the V_MA is interacting with (concentration-dependent drawdown of nutrients) is 
    V_INT<-y_farm*x_farm*t_z(t)
    # therefore nutrient change in the outflow of the farm (to depth t_z) is scaled to V_MA/V_INT of the local nutrient change (See dNH4 and dNO3 terms)
    
    
    #The effective volume that the farm interacts with over the timestep (used for phosphate calculation and also output to nutrient transport model)
    V_EFF<-y_farm*(x_farm+F_in(t))*t_z(t)
    
   
    lambda   <- (F_in(t)/x_farm) #
    
    # nutrient controls on growth, relation to biomass ####
    Q           <- ifelse(N_f>0,Q_min*(1+(N_s/N_f)),0)                                       #    Internal nutrient quota of macroalgae                                      
    B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
    g_Q         <- min(1,(Q-Q_min)/(Q-K_c))                                       # Growth limitation due to internal nutrient reserves
    # temperature-growth dynamics (from martin and marques 2002) ####
    if(SST(t)>T_O)
    {T_x<-T_max}
    else
    {T_x<-T_min}
    g_T<-exp(-2.3*((SST(t)-T_O)/(T_x-T_O))^2)
    
    # light limitation schemes####
    
    if(light_scheme==0){
      #light limits to growth including self-shading - adapted from Zollman et al 2021
      I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d(t)*h_MA+N_f*a_cs))*(1-exp(-(K_d(t)*h_MA+N_f*a_cs)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==1){
      # simple vertical light no shading
      I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d(t)*h_MA))*(1-exp(-(K_d(t)*h_MA)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==2){
      #solar angle accounted for no shading
      theta<-get_solar_angle(latitude,t)
      I_top<-PAR(t)*exp(-(K_d(t))*(z-h_MA)/sin(theta*pi/180))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d(t)*h_MA/sin(theta*pi/180)))*(1-exp(-(K_d(t)*h_MA/sin(theta*pi/180))))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==3){
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
    PO4_tot<-PO4_ext(t)*V_EFF/V_MA
    #convert N:P from mol/mol to mg/mol
    N_to_P<-N_to_P*14
    
    #capture change in N_f for harvest function
    Nf_change<-min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_f)
    
    
    #model state variables ####
    dNH4        <- lambda*(NH4_ext(t)-NH4) -((f_NH4*B)-(r_L*D)+(r_N*NH4)-(d_m*N_s))*V_MA/V_INT  # change in NH4 with time
    
    dNO3        <- lambda*(NO3_ext(t)-NO3) -((f_NO3*B)-(r_N*NH4))*V_MA/V_INT          # change in NO3 with time 
    
    dN_s        <- (f_NH4+f_NO3)*B-min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
    
    dN_f        <- min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
    
    dD          <- lambda*(D_ext(t)-D) + (d_m*N_f - r_L*D)*V_MA/V_INT                 # change in detritus with time
    
    
    dYield_farm      <- 0
  
    dYield_per_m         <- 0
    
    
    output<-list(c(dNH4, dNO3,dN_s,dN_f,dD,dYield_farm,dYield_per_m),
         Nf_change = Nf_change,
         Farm_NO3_demand = NO3_removed*V_MA,
         Farm_NH4_demand = NH4_removed*V_MA,
         V_EFF = V_EFF,
         V_INT = V_INT,
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


MA_model_const <- function(t, y, parms, data, latitude) {
  #New model version diverges from Hadley model in its spatial structure. Hadley et al model runs notionally per unit volume of water, with flow through of water in volume units. In the old version of the model implemented here this was adapted for a flow-through and vertical water flux in m/d. A major shortcoming of the old version is the resupply of nutrients via the labmda term, which cannot be greater than 1 and therefore the total input of nutrient can only be enough per timestep to restore nutrient back to ambient levels. At a timestep of 1 day this is not sufficient given the significant drawdown of nutrients by the seaweed and the model is limited then by the lambda term (the way the Hadley model is structured values of labmda>1 when F_in > farm volume increase nutrient concentration above ambient which does not make sense). In this new version of the model, advection and vertical mixing lengths give an effective volume over which nutrient concentration change is calculated, thereby allowing for a realistic throughput of nutrient and output of meaningful delta-C.
  
  #parameters expected: ####
  # Farm location:
  #  latitude (degrees)
  # Farm dimensions: 
  #  x_farm (m), 
  #  y_farm (m),
  #  z (m)
  # MA dimensions: 
  #  h_MA (m), 
  #  w_MA (m), 
  #  density_MA (fraction - 0 to 1)
  # Growth parameters: 
  #  Q_min (mg N per g dry weight), 
  #  K_c (mg N per g dry weight), 
  #  T_O (degrees celcius), 
  #  T_max (degrees celcius), 
  #  T_min (degrees celcius), 
  #  V_NH4 (mg N per g dry weight per day), 
  #  V_NO3 (mg N per g dry weight per day), 
  #  K_NH4 (mg N per m^3), 
  #  K_N03 (mg N per m^3), 
  #  a_cs  (m^2 per mg N), 
  #  I_s,  (micromol photons per m^2 per s)
  #  mu,  1/d
  #  r_L, 1/d
  #  r_N, 1/d
  #  d_m, 1/d
  # MA properties: 
  #  N_to_P, molar ratio (unitless),
  
  #time varying inputs: 
  # K_d, 1/m 
  # SST, degrees celcius
  # NO3_ext, mg N per m^3
  # NH4_ext, mg N per m^3
  # PO4_ext, mmol P per m^3 (i.e. uM)
  # D_ext,  mg N per m^3
  # t_z, m
  #F_in, m/d
  
  
  with(as.list(c(y, parms, data)), {
    
    #dimensions ####
    #farm dimensions are defined by input parms x_farm, y_farm (both 1000m i.e. grid square dimensions by default) and z_farm (depth to bottom of macroalgae)
    #Macroalgae are grown on notional 'lines' running the full width of y_farm (perpendicular to notional flow in x direction). The 'spacing' of the lines is determined by input parm density_MA. The total volume occupied by macroalgae is calculated as y_farm*x_farm*density_MA*h_MA
    V_MA<-y_farm*x_farm*density_MA*h_MA
    # The volume of water that the V_MA is interacting with (concentration-dependent drawdown of nutrients) is 
    V_INT<-y_farm*x_farm*t_z
    # therefore nutrient change in the outflow of the farm (to depth t_z) is scaled to V_MA/V_INT of the local nutrient change (See dNH4 and dNO3 terms)
    
    
    #The effective volume that the farm interacts with over the timestep (used for phosphate calculation and also output to nutrient transport model)
    V_EFF<-y_farm*(x_farm+F_in)*t_z
    
   
    lambda   <- (F_in/x_farm) #
    
    # nutrient controls on growth, relation to biomass ####
    Q           <- ifelse(N_f>0,Q_min*(1+(N_s/N_f)),0)                                       #    Internal nutrient quota of macroalgae                                      
    B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
    g_Q         <- min(1,(Q-Q_min)/(Q-K_c))                                       # Growth limitation due to internal nutrient reserves
    # temperature-growth dynamics (from martin and marques 2002) ####
    if(SST>T_O)
    {T_x<-T_max}
    else
    {T_x<-T_min}
    g_T<-exp(-2.3*((SST-T_O)/(T_x-T_O))^2)
    
    # light limitation schemes####
    
    if(light_scheme==0){
      #light limits to growth including self-shading - adapted from Zollman et al 2021
      I_top<-PAR*exp(-(K_d)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA+N_f*a_cs))*(1-exp(-(K_d*h_MA+N_f*a_cs)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==1){
      # simple vertical light no shading
      I_top<-PAR*exp(-(K_d)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA))*(1-exp(-(K_d*h_MA)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==2){
      #solar angle accounted for no shading
      theta<-get_solar_angle(latitude,t)
      I_top<-PAR*exp(-(K_d)*(z-h_MA)/sin(theta*pi/180))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA/sin(theta*pi/180)))*(1-exp(-(K_d*h_MA/sin(theta*pi/180))))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==3){
      #solar angle accounted included with shading
      theta<-get_solar_angle(latitude,t)
      sine<-sin(theta*pi/180)
      I_top<-PAR*exp(-(K_d)*(z-h_MA)/sine)
      I_av <-(I_top / ((K_d*h_MA/sine)   +   (N_f*a_cs/(sine))))*(1-exp(-(  (K_d*h_MA/sine)  +   (N_f*a_cs/(sine))   )))
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
    PO4_tot<-PO4_ext*V_EFF/V_MA
    #convert N:P from mol/mol to mg/mol
    N_to_P<-N_to_P*14
    
    #capture change in N_f for harvest function
    Nf_change<-min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_f)
    
    
    #model state variables ####
    dNH4        <- lambda*(NH4_ext-NH4) -((f_NH4*B)-(r_L*D)+(r_N*NH4)-(d_m*N_s))*V_MA/V_INT  # change in NH4 with time
    
    dNO3        <- lambda*(NO3_ext-NO3) -((f_NO3*B)-(r_N*NH4))*V_MA/V_INT          # change in NO3 with time 
    
    dN_s        <- (f_NH4+f_NO3)*B-min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
    
    dN_f        <- min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
    
    dD          <- lambda*(D_ext-D) + (d_m*N_f - r_L*D)*V_MA/V_INT                 # change in detritus with time

    output<-list(c(dNH4, dNO3, dN_s, dN_f, dD))

    return(output)
  })

  
}


MA_model_jacobian <- function(t, y, parms, data, latitude) {

  with(as.list(c(y, parms, data)), {
    
    #dimensions ####
    #farm dimensions are defined by input parms x_farm, y_farm (both 1000m i.e. grid square dimensions by default) and z_farm (depth to bottom of macroalgae)
    #Macroalgae are grown on notional 'lines' running the full width of y_farm (perpendicular to notional flow in x direction). The 'spacing' of the lines is determined by input parm density_MA. The total volume occupied by macroalgae is calculated as y_farm*x_farm*density_MA*h_MA
    V_MA<-y_farm*x_farm*density_MA*h_MA
    # The volume of water that the V_MA is interacting with (concentration-dependent drawdown of nutrients) is 
    V_INT<-y_farm*x_farm*t_z
    # therefore nutrient change in the outflow of the farm (to depth t_z) is scaled to V_MA/V_INT of the local nutrient change (See dNH4 and dNO3 terms)
    #The effective volume that the farm interacts with over the timestep (used for phosphate calculation and also output to nutrient transport model)
    V_EFF<-y_farm*(x_farm+F_in)*t_z
    
   
    lambda   <- (F_in/x_farm) #
    
    # nutrient controls on growth, relation to biomass ####
    Q           <- ifelse(N_f>0, Q_min*(1+(N_s/N_f)), 0)                                       #    Internal nutrient quota of macroalgae                                      
    B           <- N_f / Q_min                                                 # Biomass of dry macroalgae
    g_Q         <- min(1, (Q-Q_min)/(Q-K_c))                                       # Growth limitation due to internal nutrient reserves

    # temperature-growth dynamics (from martin and marques 2002) ####
    if(SST > T_O)
      {T_x <- T_max}
    else
      {T_x <- T_min}
    g_T <- exp(-2.3*((SST-T_O)/(T_x-T_O))^2)
    
    # light limitation schemes####
    
    if(light_scheme==0){
      #light limits to growth including self-shading - adapted from Zollman et al 2021
      I_top<-PAR*exp(-(K_d)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA+N_f*a_cs))*(1-exp(-(K_d*h_MA+N_f*a_cs)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==1){
      # simple vertical light no shading
      I_top <- PAR*exp(-(K_d)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA))*(1-exp(-(K_d*h_MA)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==2){
      #solar angle accounted for no shading
      theta<-get_solar_angle(latitude,t)
      I_top<-PAR*exp(-(K_d)*(z-h_MA)/sin(theta*pi/180))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA/sin(theta*pi/180)))*(1-exp(-(K_d*h_MA/sin(theta*pi/180))))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==3){
      #solar angle accounted included with shading
      theta<-get_solar_angle(latitude,t)
      sine<-sin(theta*pi/180)
      I_top<-PAR*exp(-(K_d)*(z-h_MA)/sine)
      I_av <-(I_top / ((K_d*h_MA/sine)   +   (N_f*a_cs/(sine))))*(1-exp(-(  (K_d*h_MA/sine)  +   (N_f*a_cs/(sine))   )))
      g_E <- I_av/(((I_s)+I_av)) 
    } 
    
    # uptake of nutrients from the water ####
    
    f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
    f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.

    #How much phosphate is available for growth
    PO4_tot <- PO4_ext * V_EFF/V_MA
    #convert N:P from mol/mol to mg/mol
    N_to_P <- N_to_P * 14


    ### Useful derivatives and values
    df_NH4.dNH4 = ( (V_NH4*K_NH4) / (K_NH4+NH4)^2 ) * ((Q_max-Q)/(Q_max-Q_min))
    df_NO3.dNO3 = ( (V_NO3*K_NO3) / (K_NO3+NO3)^2 ) * ((Q_max-Q)/(Q_max-Q_min))

    df_NH4.dQ = ((V_NH4*NH4)/(K_NH4+NH4)) * (1/(Q_max-Q_min))
    df_NO3.dQ = ((V_NO3*NO3)/(K_NO3+NO3)) * (1/(Q_max-Q_min))

    dQ.dNs = Q_min / N_f
    dQ.dNf = - Q_min * N_s / N_f^2

    dB.dNf = 1 / Q_min

    dg.dQ = g_T * g_E * (Q_min - K_c)/(Q - K_c)^2

    g = g_T * g_E * g_Q

    limited_P = (mu * g * N_s) > (PO4_tot * N_to_P)

    # All Jacobian values
    dJ_NH4.dNH4 = - lambda + ( - B * df_NH4.dNH4 - r_N ) * V_MA/V_INT
    dJ_NH4.dNO3 = 0
    dJ_NH4.dNs = (d_m - B * df_NH4.dQ * dQ.dNs ) * V_MA/V_INT
    dJ_NH4.dNf = ( - dB.dNf * f_NH4 - B * df_NH4.dQ * dQ.dNf ) * V_MA/V_INT
    dJ_NH4.dD = r_L * V_MA/V_INT

    dJ_NO3.dNH4 = r_N * V_MA/V_INT
    dJ_NO3.dNO3 = - lambda + ( - B * df_NO3.dNO3 )  * V_MA/V_INT
    dJ_NO3.dNs = ( - B * df_NO3.dQ * dQ.dNs ) * V_MA/V_INT
    dJ_NO3.dNf = ( - dB.dNf * f_NO3 - B * df_NO3.dQ * dQ.dNf ) * V_MA/V_INT
    dJ_NO3.dD = 0

    dJ_Ns.dNH4 = B * df_NH4.dNH4
    dJ_Ns.dNO3 = B * df_NO3.dNO3
    dJ_Ns.dNs = ( df_NH4.dQ + df_NO3.dQ ) * dQ.dNs * B - ifelse(limited_P, 0, mu * dg.dQ * dQ.dNs + mu * g) - d_m
    dJ_Ns.dNf = ( df_NH4.dQ + df_NO3.dQ ) * dQ.dNf * B - ifelse(limited_P, 0, mu * dg.dQ * dQ.dNf * N_s)
    dJ_Ns.dD = 0

    dJ_Nf.dNH4 = 0
    dJ_Nf.dNO3 = 0
    dJ_Nf.dNs = ifelse(limited_P, 0, mu * dg.dQ * dQ.dNs + mu * g)
    dJ_Nf.dNf = ifelse(limited_P, 0, mu * dg.dQ * dQ.dNf * N_s) - d_m
    dJ_Nf.dD = 0

    dJ_D.dNH4 = 0
    dJ_D.dNO3 = 0
    dJ_D.dNs = 0
    dJ_D.dNf = d_m * V_MA/V_INT
    dJ_D.dD = - lambda - r_L * V_MA/V_INT

    # !!! The result is transposed compared to its appearance here
    output = array(c( dJ_NH4.dNH4, dJ_NH4.dNO3, dJ_NH4.dNs, dJ_NH4.dNf, dJ_NH4.dD,
                      dJ_NO3.dNH4, dJ_NO3.dNO3, dJ_NO3.dNs, dJ_NO3.dNf, dJ_NO3.dD,
                      dJ_Ns.dNH4,  dJ_Ns.dNO3,  dJ_Ns.dNs,  dJ_Ns.dNf,  dJ_Ns.dD,
                      dJ_Nf.dNH4,  dJ_Nf.dNO3,  dJ_Nf.dNs,  dJ_Nf.dNf,  dJ_Nf.dD,
                      dJ_D.dNH4,   dJ_D.dNO3,   dJ_D.dNs,   dJ_D.dNf,   dJ_D.dD
                        ), dim=c(5,5))

    return(output)
  })

  
}




run_steady_state<-function(parms,...){
  
  out<-steady(y=c(NH4=1,NO3=10,N_s=100000,N_f=100000,D=100),
              time=c(0,Inf),
              func=MA_model_steady,
              parms=parms,
              method='runsteady',
              positive=TRUE)
  #make a neat array of results - if steady != 1 then steady state has failed...
  unlist(c(out[1:15],steady=attr(out,'steady')))
  
}







MA_model_steady <- function(t,y,parms,...) {
  
  
  with(as.list(c(y, parms)), {
    
    V_MA<-y_farm*x_farm*density_MA*h_MA
    V_INT<-y_farm*(x_farm)*(t_z)
    V_EFF<-y_farm*(x_farm+F_in)*(t_z)
    
    # therefore nutrient change in the outflow of the farm (to depth z+t_z) is scaled to V_MA/V_EFF of the local nutrient change (See dNH4 and dNO3 terms)
    
    lambda   <- F_in/x_farm #
    
    # nutrient controls on growth, relation to biomass ####
    Q           <- ifelse(N_f>0,Q_min*(1+(N_s/N_f)),0)                                       #    Internal nutrient quota of macroalgae                                      
    B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
    g_Q         <- min(1,(Q-Q_min)/(Q-K_c)) #min(1,((Q-Q_min)/(Q-K_c)))                                      # Growth limitation due to internal nutrient reserves
    # temperature-growth dynamics (from martin and marques 2002) ####
    if(SST>T_O)
    {T_x<-T_max}
    else
    {T_x<-T_min}
    g_T<-exp(-2.3*((SST-T_O)/(T_x-T_O))^2)
    
    # light limitation schemes####
    
    if(light_scheme==0){
      #light limits to growth including self-shading - adapted from Zollman et al 2021
      I_top<-PAR*exp(-(K_d)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA+N_f*a_cs))*(1-exp(-(K_d*h_MA+N_f*a_cs)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==1){
      # simple vertical light no shading
      I_top<-PAR*exp(-(K_d)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA))*(1-exp(-(K_d*h_MA)))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==2){
      #solar angle accounted for no shading
      theta<-get_solar_angle(latitude,doy=180)
      I_top<-PAR*exp(-(K_d)*(z-h_MA)/sin(theta*pi/180))                                    # calculate incident irradiance at the top of the farm
      I_av <-(I_top/(K_d*h_MA/sin(theta*pi/180)))*(1-exp(-(K_d*h_MA/sin(theta*pi/180))))          # calclulate the average irradiance throughout the height of the farm
      g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
    } else if (light_scheme==3){
      #solar angle accounted included with shading
      theta<-get_solar_angle(latitude,doy=180)
      sine<-sin(theta*pi/180)
      I_top<-PAR*exp(-(K_d)*(z-h_MA)/sine)
      I_av <-(I_top / ((K_d*h_MA/sine)   +   (N_f*a_cs/(sine))))*(1-exp(-(  (K_d*h_MA/sine)  +   (N_f*a_cs/(sine))   )))
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
    PO4_tot<-PO4_ext*V_EFF/V_MA
    #convert N:P from mol/mol to mg/mol
    N_to_P<-N_to_P*14
    
    
    
    
    #model state variables ####
    dNH4        <- lambda*(NH4_ext-NH4) -((f_NH4*B)-(r_L*D)+(r_N*NH4)-(d_m*N_s))*V_MA/V_INT  # change in NH4 with time
    
    dNO3        <- lambda*(NO3_ext-NO3) -((f_NO3*B)-(r_N*NH4))*V_MA/V_INT           # change in NO3 with time 
    
    dN_s        <- (f_NH4+f_NO3)*B-min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
    
    dN_f        <- min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
    
    dD          <- lambda*(D_ext-D) + (d_m*N_f - r_L*D)*V_MA/V_INT               # change in detritus with time
    
    
    
    
    
    output<-list(c(dNH4, dNO3,dN_s,dN_f,dD),
                 Farm_NO3_demand = NO3_removed*V_MA,
                 Farm_NH4_demand = NH4_removed*V_MA,
                 V_EFF = V_EFF,
                 V_INT=V_INT,
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






