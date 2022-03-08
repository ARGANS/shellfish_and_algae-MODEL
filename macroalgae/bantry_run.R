### test run with bantry bay data - MJ 1-/02/2022

source('run_MA.R')
source('calc_par.R')
library(atmos)



bantry<-read.csv('bantry_data/bantry_MLDaveraged.csv',sep = ';')


#the new data will already have unit conversion done..
#bantry<-read.csv('bantry_data/bantry_3m.csv',sep = ';')

#quick hack substitute missing values
bantry$Ammonium[is.na(bantry$Ammonium)]<-0.4
bantry$Nitrate[is.na(bantry$Nitrate)]<-8



convert_umolN_to_mgm3<-function(data_in){
  rmm_N<-14 #gmol
  mmol_m3<-data_in #(1mmol m^(-3) = 1umol l-1)
  mg_m3<-data_in*14
  mg_m3
}

convert_photonflux<-function(data_in){
  #converts from molm-2d-1 to umolm-2s-1
  data_in*1e6/(3600*24)
}

bantry$Ammonium<-convert_umolN_to_mgm3(bantry$Ammonium)
bantry$Nitrate<-convert_umolN_to_mgm3(bantry$Nitrate)
bantry$F_in<-sqrt(bantry$northward_Water_current^2 + bantry$eastward_Water_current^2)*3600*24 #take 'hypotenuse of N and E flow and convert from m/s to m/day
bantry$F_in[is.na(bantry$F_in)]<-2500

bantry$par<-convert_photonflux(bantry$par) 

##used for par substitution documentation only
obspar<-bantry$par
##
par_GC<-read.csv('bantry_data/bantry_par_theoretical.csv')

bantry$par[1:365]<-substitute_par(bantry$par[1:365],par_GC$par_GC)
bantry$par[366:396]<-bantry$par[1:31]

#make a 2year input data
byear<-bantry[1:365,]
b2y<-rbind(byear,byear)


bantry_input_data<-data.frame(
  time=1:730,
  PAR=b2y$par,
  SST=b2y$Temperature,
  NH4_in=b2y$Ammonium,
  NO3_in=b2y$Nitrate,
  PO4_in=1000,
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
  light_scheme=4,
  latitude=51.6

)





bantry_in<-setup_run_input(input=bantry_input_data)
bantry_parms<-setup_run_parms(parms=c(parms_bantry_alaria,parms_bantry))
bantry_springharvest_parms<-setup_run_parms(parms=c(parms_bantry_alaria,parms_bantry,harvest_winter_growth_run,harvest_method=1))
bantry_continuousharvest_parms<-setup_run_parms(parms=c(parms_bantry_alaria,parms_bantry,harvest_CCA_run,harvest_method=2))


bantry_run<-run_MA_model(input=bantry_in,
                         parameters=bantry_parms,
                         y0=c(NH4=bantry$Ammonium[1],NO3=bantry$Nitrate[1],N_s=100,N_f=100,D=0,Yield_farm=0, Yield_per_m=0)
                         )


bantry_spring_harvest_run<-run_MA_model(input=bantry_in,
                         parameters=bantry_springharvest_parms,
                         y0=c(NH4=bantry$Ammonium[1],NO3=bantry$Nitrate[1],N_s=0,N_f=0,D=0,Yield_farm=0, Yield_per_m=0)
                         )

bantry_CCA_run<-run_MA_model(input=bantry_in,
                                        parameters=bantry_continuousharvest_parms,
                                        y0=c(NH4=bantry$Ammonium[1],NO3=bantry$Nitrate[1],N_s=0,N_f=0,D=0,Yield_farm=0, Yield_per_m=0)
)


#names(bantry_run)[names(bantry_run)=='time'][2]<-'time2'
library(reshape2)
library(ggplot2)
library(patchwork)
plot_bantry_results<-function(x.,harvest=FALSE){
  
  
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
