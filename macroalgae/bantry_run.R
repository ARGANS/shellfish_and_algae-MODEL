### test run with bantry bay data - MJ 1-/02/2022

bantry<-read.csv('bantry_data/bantry_MLDaveraged_2020.csv',sep = ';')

source('run_MA.R')


#assume nutrient data in umol - to check

convert_umolN_to_mgm3<-function(data_in){
  rmm_N<-14 #gmol
  mmol_m3<-data_in #(1mmol m^(-3) = 1umol l-1)
  mg_m3<-data_in*rmm_N
  mg_m3
}

convert_photonflux<-function(data_in){
  #converts from molm-2d-1 to umolm-2s-1
  data_in*1e6/(3600*24)
}



bantry$Ammonium<-convert_umolN_to_mgm3(bantry$Ammonium)
bantry$Nitrate<-convert_umolN_to_mgm3(bantry$Nitrate)
bantry$F_in<-sqrt(bantry$northward_Water_current^2 + bantry$eastward_Water_current^2)*3600*24 #take 'hypotenuse of N and E flow and convert from m/s to m/day
bantry$par<-convert_photonflux(bantry$par) 


##### idealised PAR for substitution
bantryGCPAR<-data.frame(time=1:365)
latitude_bantry<-51.9
longitude_bantry<-9.8


####### wrapper function to run Greg Carder model for PAR and get the result we need out
library(atmos)

run_GC.f<-function(jday,lat,long,hr){
   GCOutput<-GreggCarder.f(jday = jday,rlon = long,rlat=lat,hr = hr)

  
  ###### convert GC model output to total par in umol m-2 s-1
  Convert_GC_ED<-function(GCOutput){
    ED<-GCOutput$Ed
    LAM<-GCOutput[['lam']]
    #these terms from GregCarder.f
    h <-6.6256e-34
    c <- 299800000
    hc<-1/(h*c)
    sumpar<-1e-09* hc/6.023e+17 * sum(ED*LAM)
    sumpar
    
  }
  Convert_GC_ED(GCOutput)
}


daily_GC.f<-function(lat,long,jday){
  #function to calculate day-averaged photon flux
  #longitude is not needed (only relevant for instantaneous flux at given time)
  #we only need to run this for each latitude band in the model. 
  timeseq<-0:23
  sum(sapply(0:23,run_GC.f,long=long,lat=lat,jday=jday),na.rm=T)/24
  
}

bantryGCPAR$PAR<-sapply(bantryGCPAR$time,daily_GC.f,lat=latitude_bantry,long=longitude_bantry)


input_data<-data.frame(
  time=1:345,
  PAR=bantry$par,
  SST=bantry$Temperature,
  NH4_in=bantry$Ammonium,
  NO3_in=bantry$Nitrate,
  PO4_in=1000,
  K_d=0.11,
  F_in=bantry$F_in,
  theta = setup_solar_angle(latitude=51.7,start_day=25,ndays=345)
)



#input_data_OctStart<-rbind(input_data[330:345],)

parms_bantry_alaria <- c(
  mu      = 0.4,    # maximum growth rate           / 1/d
  V_NH4   = 60,      # max. ammonium uptake rate     / mg(N)g-1(dw)d-1
  V_NO3   = 200,      # max. nitrate uptake rate      / mg(N)g-1(dw)d-1
  K_NH4   = 350,     # Half saturation constant NH4  / mg(N)m-3
  K_NO3   = 500,     #   "      "          "    NO3  / mg(N)m-3
  Q_max   = 50,      # max. internal nitrogen        / mg(N) g-1 (dw)
  Q_min   = 14,      # min.    "        "            / mg(N) g-1 (dw)
  N_to_P  = 12,      #N:P ratio of seaweed biomass
  K_c     = 7,       # Half growth constant          / mg(N) g-1 (dw)
  T_O     = 8,      # optimum growth temperature    / oC
  T_min     = 1,       # min temperature for growth  / oC
  T_max    = 19,     # max temperature for growth     / oC
  T_r     = 1,
  I_s     = 150,     # saturation irradiance         / umol photons m-2 s-1
  a_cs    = 0.00036, # nitrogen-specific shading     / m2 mg-1 (N)
  d_m     = 0.003,   # mortality rate                / 1/d
  h_MA    = 1.2,     # height of seaweed             / m
  w_MA    = 0.2,     # width of seaweed e.g. on rope /m
  r_L     = 0.2,     # remineralisation rate         / 1/d
  r_N     = 0.1     # nitrification rate            / 1/d
)


bantry_run<-run_model(c(default_parms_farm,default_parms_ulva,default_parms_run),default_input,input=input_data,parms=parms_bantry_alaria,y0=c(c(NH4=bantry$Ammonium[1],NO3=bantry$Nitrate[1],N_s=100,N_f=100,D=0,Yield=0)))
names(bantry_run)[names(bantry_run)=='time'][2]<-'time2'

library(ggplot2)

ggplot(bantry_run) +
  aes(x = time, y = PAR) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal()

default_input %>%
  filter(time >= 0L & time <= 400L) %>%
  ggplot() +
  aes(x = time, y = PAR) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal()

plot(bantry_run$time,bantry_run$B,xlab='day of year',ylab='Dry Weight of Macroalgae /g')
#plot(bantry_run$time,bantry_run$Yield/parms_bantry_alaria[['Q_min']],xlab='day of year',ylab='Dry Weight of Macroalgae cumulative yield /g')
