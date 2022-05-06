#run_SF.R - controls for running model manually from R. Reference runs for documentation

source('shellfish_model.R')




# dummy data (in lieu of environmental forcings) ######## -------------------------

# construct dummy input data

### temperature  # celcius
SST_mean <- 10
SST_magn <- 10

POC_mean <- 500
POC_magn <- 300

CHL_mean <- 6
CHL_magn  <-5

PHYC_mean<-95
PHYC_magn<-75



run_length <- time <- 1:(2*365)



default_input<- data.frame(
  time   = time,
  
  SST    = SST_mean + (SST_magn*sin(2*pi*(time+180)/365)),
  CHL    = CHL_mean + (CHL_magn*sin(2*pi*(time+180)/365)),
  PHYC  = PHYC_mean + (PHYC_magn*sin(2*pi*(time+180)/365)),
  POC    = POC_mean + (POC_magn*sin(2*pi*(time+180)/365)),
  
  F_in   = 1000,
  t_z    = 10
)

##load in all parameters from json file
library(rjson)
allparams<-fromJSON(file='shellfish_model_parameters.json')

parms_M_edulis<-unlist(allparams$species$defaults$M_edulis$parameters)
parms_C_gigas<-unlist(allparams$species$defaults$C_gigas$parameters)
parms_P_maximus<-unlist(allparams$species$defaults$P_maximus$parameters)


default_parms_SF_farm<-c(unlist(allparams$farm$defaults$default$parameters),unlist(allparams$farm$defaults$default$options))
default_run_parms<-c(unlist(allparams$run$defaults$default$parameters),c(allparams$run$defaults$default$options))

## make full default parameter set --------
default_parms<-c(default_run_parms,default_parms_SF_farm,parms_M_edulis)

setup_run_parms<-function(default_parms.=default_parms, parms){
  #when calling rum_SF_model from R, use this function to overwrite default parameters as needed
  run_parms<-replace(default_parms.,names(parms),parms)
}

setup_run_input<-function(default_input.=default_input,input){
  df<-default_input.[1:nrow(input),]
  df[,colnames(input)]<-input
  df
}

# Idealised reference runs
#-------------------------------

# These are run using the reference_run() function below by applying different input_data as defined here for the 4 reference runs. 

#define 'low' and 'high' values for each parameter
low_chl_mean<-3
low_chl_magn<-2.9
high_chl_mean<-20
high_chl_magn<-13
low_POC_mean<-100
low_POC_magn<-50
high_POC_mean<-500
high_POC_magn<-250
low_flow<-10
low_flow_tz<-5
high_flow<-1000
high_flow_tz<-10

##1 Low chlorophyll, low POC, low flow

lclf<-c(
  CHL_mean=low_chl_mean,
  CHL_magn=low_chl_magn,
  POC_mean=low_POC_mean,
  POC_magn=low_POC_magn,
  F_in=low_flow,
  t_z=low_flow_tz
)


lchf<-c(
  CHL_mean=low_chl_mean,
  CHL_magn=low_chl_magn,
  POC_mean=low_POC_mean,
  POC_magn=low_POC_magn,
  F_in=high_flow,
  t_z=high_flow_tz
)


hchf<-c(
  CHL_mean=high_chl_mean,
  CHL_magn=high_chl_magn,
  POC_mean=high_POC_mean,
  POC_magn=high_POC_magn,
  F_in=high_flow,
  t_z=high_flow_tz
)

hclf<-c(
  CHL_mean=high_chl_mean,
  CHL_magn=high_chl_magn,
  POC_mean=high_POC_mean,
  POC_magn=high_POC_magn,
  F_in=low_flow,
  t_z=low_flow_tz
)

reference_run <- function(input_data,nondefault_parms,species='Med'){
  #to run with other species, 'Med', 'Cgig', 'Pmax'
  with(as.list(input_data), {
    time<-1:(3*365)
    input_frame<- data.frame(
      time   = time,
      SST    = 13 + (4*sin(2*pi*(time+180)/365)),
      CHL = CHL_mean + (CHL_magn*sin(2*pi*(time+91)/365)),
      POC = POC_mean + (POC_magn*sin(2*pi*(time+91)/365)),
      F_in   = F_in,
      t_z    = t_z
    )
    
    if(species=='Cgig'){
      nondefault_parms=c(nondefault_parms,parms_C_gigas)
    } else if(species=='Pmax'){
      nondefault_parms=c(nondefault_parms,parms_P_maximus)
    }
    parms.=setup_run_parms(parms=nondefault_parms)
    y0=c(SHE=100,STE=1000,spawnday=365,sd2=365,CHL_farm=input_frame$CHL[1],POC_farm=input_frame$POC[1],POP=parms.$stocking_density,PHYC_farm=0)
    run_SF_model(input=input_frame,parameters = parms.,y0=y0,output='df')
  })
  
}

HL<-reference_run(hclf,c(PHYC_data=FALSE))
LL<-reference_run(lclf,c(PHYC_data=FALSE))
LH<-reference_run(lchf,c(PHYC_data=FALSE))
HH<-reference_run(hchf,c(PHYC_data=FALSE))


