### test run with bantry bay data - MJ 1-/02/2022

source('run_MA.R')
source('calc_par.R')
library(atmos)

#latitude at bantry
bantrylat<-51.6

#bantry_old<-read.csv('bantry_data/bantry_MLDaveraged.csv',sep = ';')


#the new data will already have unit conversion done..
bantry<-read.csv('bantry_data/bantry_3m.csv',sep = ';')

#quick hack substitute missing values
bantry$Ammonium[is.na(bantry$Ammonium)]<-4
bantry$Nitrate[is.na(bantry$Nitrate)]<-80
bantry$Temperature[is.na(bantry$Temperature)]<-8


#convert_umolN_to_mgm3<-function(data_in){
#   rmm_N<-14 #gmol
#   mmol_m3<-data_in #(1mmol m^(-3) = 1umol l-1)
#   mg_m3<-data_in*14
#   mg_m3
# }
# 
# convert_photonflux<-function(data_in){
#   #converts from molm-2d-1 to umolm-2s-1
#   data_in*1e6/(3600*24)
# }
# 
# bantry$Ammonium<-convert_umolN_to_mgm3(bantry$Ammonium)
# bantry$Nitrate<-convert_umolN_to_mgm3(bantry$Nitrate)
# bantry$F_in<-sqrt(bantry$northward_Water_current^2 + bantry$eastward_Water_current^2)*3600*24 #take 'hypotenuse of N and E flow and convert from m/s to m/day
bantry$F_in<-bantry$current_intensity
bantry$F_in[is.na(bantry$F_in)]<-2500

#bantry$par<-convert_photonflux(bantry$par) 

##used for par substitution documentation only
obspar<-bantry$par
##
#par_GC<-read.csv('bantry_data/bantry_par_theoretical.csv')
#PARsub<-par_GC$par_GC


PARsub<-get_simulated_PAR_by_latitdue(bantrylat)




bantry$par[1:365]<-substitute_par(bantry$par[1:365],PARsub)
bantry$par[366:396]<-bantry$par[1:31]

#make a 2year input data
byear<-bantry[1:365,]
b2y<-rbind(byear,byear)


bantry_input_data<-data.frame(
  time=1:730,
  PAR=b2y$par,
  SST=b2y$Temperature,
  NH4_ext=b2y$Ammonium,
  NO3_ext=b2y$Nitrate,
  PO4_ext=1000,
  K_d=0.1,
  F_in=b2y$F_in,
  t_z = 10
)



parms_bantry_alaria <- c(
  mu      = 0.12,    # maximum growth rate           / 1/d
  V_NH4   = 100,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
  V_NO3   = 200,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
  K_NH4   = 1,     # Half saturation constant NH4  / mg(N)m-3
  K_NO3   = 2,     #   "      "          "    NO3  / mg(N)m-3
  Q_max   = 22,      # max. internal nitrogen        / mg(N) g-1 (dw)
  Q_min   = 10,      # min.    "        "            / mg(N) g-1 (dw)
  N_to_P  = 12,      #N:P ratio of seaweed biomass
  K_c     = 8,       # Half growth constant          / mg(N) g-1 (dw)
  T_O     = 12.5,      # optimum growth temperature    / oC
  T_min     = 0,       # min temperature for growth  / oC
  T_max    = 20,     # max temperature for growth     / oC
  T_r     = 1,
  I_s     =90,     # saturation irradiance         / umol photons m-2 s-1
  a_cs    = 0.00036, # nitrogen-specific shading     / m2 mg-1 (N)
  d_m     = 0.002,   # mortality rate                / 1/d
  h_MA    = 2,     # height of seaweed             / m
  w_MA    = 0.3,     # width of seaweed e.g. on rope /m
  r_L     = 0.10,     # remineralisation rate         / 1/d
  r_N     = 0.1     # nitrification rate            / 1/d
)

parms_bantry<-c(
  light_scheme=3,
  latitude=bantrylat

)





bantry_in<-setup_run_input(input=bantry_input_data)
bantry_parms<-setup_run_parms(parms=c(parms_bantry_alaria,parms_bantry))
  parmsjson<-toJSON(as.list(bantry_parms))
  write(parmsjson,file = 'Bantry_run_model_parms.json')
bantry_springharvest_parms<-setup_run_parms(parms=c(parms_bantry_alaria,parms_bantry,harvest_winter_growth_run,harvest_method=1))
bantry_continuousharvest_parms<-setup_run_parms(parms=c(parms_bantry_alaria,parms_bantry,harvest_CCA_run,harvest_method=2))


bantry_run<-run_MA_model(input=bantry_in,
                         parameters=bantry_parms,
                         y0=c(NH4=bantry_in$NH4_ext[1],NO3=bantry_in$NO3_ext[1],N_s=100,N_f=100,D=0,Yield_farm=0, Yield_per_m=0)
                         )


bantry_spring_harvest_run<-run_MA_model(input=bantry_in,
                         parameters=bantry_springharvest_parms,
                         y0=c(NH4=bantry_in$NH4_ext[1],NO3=bantry_in$NO3_ext[1],N_s=0,N_f=0,D=0,Yield_farm=0, Yield_per_m=0)
                         )

bantry_CCA_run<-run_MA_model(input=bantry_in,
                                        parameters=bantry_continuousharvest_parms,
                                        y0=c(NH4=bantry_in$NH4_ext[1],NO3=bantry_in$NO3_ext[1],N_s=0,N_f=0,D=0,Yield_farm=0, Yield_per_m=0)
)

#write.csv(bantry_run,'bantry_data/bantry_run.csv')
#write.csv(bantry_spring_harvest_run,'bantry_data/bantry_spring_harvest_run.csv')
#write.csv(bantry_CCA_run,'bantry_data/bantry_CCA_run.csv')

#names(bantry_run)[names(bantry_run)=='time'][2]<-'time2'
library(reshape2)
library(ggplot2)
library(patchwork)
plot_bantry_results<-function(x.,harvest=FALSE,indata=bantry_in){
  
  x.<-cbind(indata,x.)
  #x. is output of run_MA_model
  #melt it to give overlayable plots as needed
  TS<-5 #text size control
  
  ammonium<-ggplot(melt(x.[,c(1,4,12)],id='time')) +
    geom_point(shape = "circle", size = 0.7,aes(x = time, y = value,colour=variable)) +
    scale_color_hue(direction = 1) + ylab('ammonium concentration / mg N m^-3')+
    theme_minimal()+
    theme(text = element_text(size = TS)) + theme(legend.position="bottom") 
    
  nitrate<-ggplot(melt(x.[,c(1,5,13)],id='time')) +
    geom_point(shape = "circle", size = 0.7,aes(x = time, y = value,colour=variable)) +
    scale_color_hue(direction = 1) + ylab('nitrate concentration / mg N m^-3')+
    theme_minimal()+
    theme(text = element_text(size = TS)) + theme(legend.position="bottom") 
    
  growth<-ggplot(melt(x.[,c(1,27:30)],id='time')) +
    geom_point(shape = "circle", size = 0.7,aes(x = time, y = value,colour=variable)) +
    scale_color_hue(direction = 1) +  ylab('growth parameters \n / unitless except mu_g_eqt in per day')+
    theme_minimal()+
    theme(text = element_text(size = TS)) + theme(legend.position="bottom") 
    
  Ndemand<-ggplot(melt(x.[,c(1,20,21)],id='time')) +
    geom_point(shape = "circle", size = 0.7,aes(x = time, y = value/1e6,colour=variable)) +
    scale_color_hue(direction = 1) +
    ylab('N demand / tonnes')+ theme(legend.position="bottom")+
    theme_minimal()+
    theme(text = element_text(size = TS)) + theme(legend.position="bottom") 
    
  #convert units
  x.$B_line<-x.$B_line/1000  #kg
  x.$V_EFF<-x.$V_EFF/1e9 #0.1km3
  x.$Yield_farm<-x.$Yield_farm/1e9  #ktonnes
  x.$Yield_per_m<-x.$Yield_per_m/1000 #kg
  x.$F_in<-x.$F_in/1000
  Biomass<-ggplot(melt(x.[,c(1,26)],id='time')) +
    geom_point(shape = "circle", size = 0.7,aes(x = time, y = value),colour='firebrick') +
    scale_color_hue(direction = 1) +
    theme_minimal()+ylab('Biomass / (kg DW)/m')+ theme(text = element_text(size = TS))
  VEFF<-ggplot(melt(x.[,c(1,8)],id='time')) +
    geom_point(shape = "circle", size = 0.7,aes(x = time, y = value,colour=variable)) +
    scale_color_hue(direction = 1) +
    theme_minimal()+ylab('F_in / km')+ theme(text = element_text(size = TS))+
    theme(legend.position="bottom") 
  Yield<-ggplot(melt(x.[,c(1,17,18)],id='time')) +
    geom_point(shape = "circle", size = 0.7,aes(x = time, y = value,colour=variable)) +
    scale_color_hue(direction = 1) +ylab('Farm yield / ktonne \n Yield per metre / kg') +
    theme_minimal()+
    theme(text = element_text(size = TS)) + theme(legend.position="bottom") 
    
  
  if(harvest){
    VEFF+ammonium+nitrate+growth+Biomass+Yield+plot_layout(ncol = 2)
  } else {
    VEFF+ammonium+nitrate+growth+Ndemand+Biomass+ plot_layout(ncol = 2)
  }
  
}


#ggplot(bantry_run) +
#  aes(x = time, y = PAR) +
#  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
#  theme_minimal()

#default_input %>%
#  filter(time >= 0L & time <= 400L) %>%
#  ggplot() +
#  aes(x = time, y = PAR) +
#  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
#  theme_minimal()


#plot(bantry_run$time,bantry_run$B,xlab='day of year',ylab='Dry Weight of Macroalgae /g')





###### experimenting with rootsolve on monthly mean values. 




#add a date column
library(dplyr)
library(lubridate)
b2y$date<-as.Date("2020-12-31")+(1:nrow(b2y))


b2ym<- b2y %>% group_by(month=floor_date(date,"month")) %>% 
  summarize(Ammonium=mean(Ammonium,na.rm=TRUE),
            Nitrate=mean(Nitrate,na.rm=TRUE),
            Temperature=mean(Temperature,na.rm=TRUE),
            par=mean(par,na.rm=TRUE),
            F_in=mean(F_in,na.rm=TRUE)
            )

bym<-as.data.frame(b2ym)[1:12,]  

bym$month<-c(1:12)





bantry_monthly_input_data<-data.frame(
  time=bym$month,
  PAR=bym$par,
  SST=bym$Temperature,
  NH4_ext=bym$Ammonium,
  NO3_ext=bym$Nitrate,
  PO4_ext=1000,
  K_d=0.1,
  F_in=bym$F_in,
  t_z = 10,
  D_ext=0
)



run_monthly_steady_state<-function(month, parms,...){
  #get a steady state solution for any month (1-12) in the input data 
  x<-run_steady_state(parms=c(parms,unlist(bantry_monthly_input_data[month,2:10])))
  x
}



#MA_model_steady <- function(t,y,parms,...) {
#  
#   
#   with(as.list(c(y, parms)), {
#     
#      V_MA<-y_farm*x_farm*density_MA*h_MA
#      V_INT<-y_farm*(x_farm)*(t_z)
#      V_EFF<-y_farm*(x_farm+F_in)*(t_z)
#      
#     # therefore nutrient change in the outflow of the farm (to depth z+t_z) is scaled to V_MA/V_EFF of the local nutrient change (See dNH4 and dNO3 terms)
#     
#     lambda   <- F_in/x_farm #
#     
#     # nutrient controls on growth, relation to biomass ####
#     Q           <- ifelse(N_f>0,Q_min*(1+(N_s/N_f)),0)                                       #    Internal nutrient quota of macroalgae                                      
#     B           <- N_f/Q_min                                                 # Biomass of dry macroalgae
#     g_Q         <- min(1,(Q-Q_min)/(Q-K_c)) #min(1,((Q-Q_min)/(Q-K_c)))                                      # Growth limitation due to internal nutrient reserves
#     # temperature-growth dynamics (from martin and marques 2002) ####
#     if(SST>T_O)
#     {T_x<-T_max}
#     else
#     {T_x<-T_min}
#     g_T<-exp(-2.3*((SST-T_O)/(T_x-T_O))^2)
#     
#     # light limitation schemes####
#     
#     if(light_scheme==0){
#       #light limits to growth including self-shading - adapted from Zollman et al 2021
#       I_top<-PAR*exp(-(K_d)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
#       I_av <-(I_top/(K_d*h_MA+N_f*a_cs))*(1-exp(-(K_d*h_MA+N_f*a_cs)))          # calclulate the average irradiance throughout the height of the farm
#       g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
#     } else if (light_scheme==1){
#       # simple vertical light no shading
#       I_top<-PAR*exp(-(K_d)*(z-h_MA))                                    # calculate incident irradiance at the top of the farm
#       I_av <-(I_top/(K_d*h_MA))*(1-exp(-(K_d*h_MA)))          # calclulate the average irradiance throughout the height of the farm
#       g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
#     } else if (light_scheme==2){
#       #solar angle accounted for no shading
#       theta<-get_solar_angle(latitude,doy=180)
#       I_top<-PAR*exp(-(K_d)*(z-h_MA)/sin(theta*pi/180))                                    # calculate incident irradiance at the top of the farm
#       I_av <-(I_top/(K_d*h_MA/sin(theta*pi/180)))*(1-exp(-(K_d*h_MA/sin(theta*pi/180))))          # calclulate the average irradiance throughout the height of the farm
#       g_E <- I_av/(((I_s)+I_av))                                          # light limitation scaling function
#     } else if (light_scheme==3){
#       #solar angle accounted included with shading
#       theta<-get_solar_angle(latitude,doy=180)
#       sine<-sin(theta*pi/180)
#       I_top<-PAR*exp(-(K_d)*(z-h_MA)/sine)
#       I_av <-(I_top / ((K_d*h_MA/sine)   +   (N_f*a_cs/(sine))))*(1-exp(-(  (K_d*h_MA/sine)  +   (N_f*a_cs/(sine))   )))
#       g_E <- I_av/(((I_s)+I_av)) 
#     } 
#     
#     
#     # growth function (uptake of nutrient from internal reserves) ####
#     
#     mu_g_EQT    <- mu*g_E*g_Q*g_T                                            # Growth function for macroalgae
#     
#     
#     # uptake of nutrients from the water ####
#     
#     f_NH4       <- ((V_NH4*NH4)/(K_NH4+NH4))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NH4
#     f_NO3       <- ((V_NO3*NO3)/(K_NO3+NO3))*((Q_max-Q)/(Q_max-Q_min))          # uptake rate of NO3.
#     
#     # amount of nutrient removed per m3 
#     NO3_removed <- (f_NO3*B)-(r_N*NH4)
#     NH4_removed <- (f_NH4*B)-(r_L*D)+(r_N*NH4)-(d_m*N_s)
#     
#     #How much phosphate is available for growth
#     PO4_tot<-PO4_in*V_EFF/V_MA
#     #convert N:P from mol/mol to mg/mol
#     N_to_P<-N_to_P*14
#     
#     
#     
#     
#     #model state variables ####
#     dNH4        <- lambda*(NH4_in-NH4) -((f_NH4*B)-(r_L*D)+(r_N*NH4)-(d_m*N_s))*V_MA/V_INT  # change in NH4 with time
#     
#     dNO3        <- lambda*(NO3_in-NO3) -((f_NO3*B)-(r_N*NH4))*V_MA/V_INT           # change in NO3 with time 
#     
#     dN_s        <- (f_NH4+f_NO3)*B-min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_s)                          # change in internal nitrogen store - eq 5 in Hadley
#     
#     dN_f        <- min(mu_g_EQT*N_f,PO4_tot*N_to_P)-(d_m*N_f)                                    # change in fixed nitrogen (i.e. biomass nitrogen) - eq 6 in Hadley
#     
#     dD          <- lambda*(D_in-D) + (d_m*N_f - r_L*D)*V_MA/V_INT               # change in detritus with time
#     
#     
#     
#     
#     
#     output<-list(c(dNH4, dNO3,dN_s,dN_f,dD),
#                  Farm_NO3_demand = NO3_removed*V_MA,
#                  Farm_NH4_demand = NH4_removed*V_MA,
#                  V_EFF = V_EFF,
#                  V_INT=V_INT,
#                  Q         = Q,
#                  B_vol  = B,
#                  B_line = B*h_MA*w_MA,
#                  g_Q       = g_Q,
#                  g_T       = g_T,
#                  g_E       = g_E,
#                  mu_g_EQT  = mu_g_EQT,
#                  f_NH4     = f_NH4,
#                  f_NO3     = f_NO3
#     )
#   }) 
# }

#example uses
#steady(y=c(NH4=1,NO3=10,N_s=10000,N_f=10000,D=100),time=c(0,Inf),func=MA_model_steady,parms=c(default_parms,unlist(bantry_monthly_input_data[1,2:10])),method='runsteady',positive=TRUE)
#run_steady_state(parms=c(default_parms,unlist(bantry_monthly_input_data[1,2:10])))

monthly_bantry<-cbind(data.frame(month=factor(month.abb,levels=month.abb),as.data.frame(t(sapply(seq(1:12),run_monthly_steady_state,parms=default_parms)))))

plot_monthly_data<-function(x){
  ggplot(x,aes(x=month,y=B_line/1000))+geom_col()+ylab('Biomass (g DW/metre)')
}

plot_monthly_data(monthly_bantry)




### run prognostic model with daily values averaged to month

days<-seq(from=dmy('01/01/2021'),to=dmy('31/12/2021'),length.out=365)

month(days)

monthlyaveraged_daily<-bantry_monthly_input_data[1,][-1,]
for(daymonth in month(days)){
  monthlyaveraged_daily<-rbind(monthlyaveraged_daily,bantry_monthly_input_data[daymonth,])
}

MAD<-rbind(monthlyaveraged_daily,monthlyaveraged_daily)

MAD$time<-1:730

bantry_MAD_in<-setup_run_input(input=MAD)

bantry_MAD_run<-run_MA_model(input=bantry_MAD_in,
                         parameters=bantry_parms,
                         y0=c(NH4=bantry_in$NH4_ext[1],NO3=bantry_in$NO3_ext[1],N_s=100,N_f=100,D=0,Yield_farm=0, Yield_per_m=0)
)


# try with no flow

MADF<-MAD
MADF$F_in<-0
bantry_MADF_in<-setup_run_input(input=MADF)

bantry_MADF_run<-run_MA_model(input=bantry_MADF_in,
                             parameters=bantry_parms,
                             y0=c(NH4=bantry_in$NH4_ext[1],NO3=bantry_in$NO3_ext[1],N_s=100,N_f=100,D=0,Yield_farm=0, Yield_per_m=0)
)
