### Shellfish model for EMFF shellfish and algae project
## Based  on the shellfish component of Dabrowski et al 2013 http://dx.doi.org/10.1016/j.seares.2012.10.012Contents 

## Full details of scientific rationale can be found in associated ATBD
# Authors
# M. Johnson (Bantry Marine Research Station, Ireland) mjohnson@bmrs.ie,


# REQUIRES R>4.1!




#### load libraries ######

library('deSolve')
library('rootSolve')
library('tidyverse')

boundary_forcings<-function(input_data){
  # This function takes input environmental data for an individual grid square
  # as a data frame of daily average values and generates the necessary functions
  # to return interpolated values at any time point within the range
  # (i.e. at whatever temporal resolution is required  by the ODE solver).
  # The following columns are expected in the data frame. 
  
  # time (in days, e.g. 1:365),
  # SST (Sea surface temperature in celcius),
  # Chlorophyll concentration (in g m-3 [3])
  # F_in (net horizontal flow rate into farm cross section (m/d))
  # t_z (vertical mixing depth (m/d))
  
  output = list()
  for (param in names(input_data)) {
    output[[param]] = approxfun(x=input_data$time, y=input_data[,param], method='linear', rule=2)
  }
  
  return(output)

}








## =============================================================================
# event control functions - when to plant out and harvest the shellfish and trigger spawning events
## =============================================================================

spawn_timed_rootfunc<-function(t,y,parms,...){
  with(as.list(c(y, parms)),{
    SHL<-a_SHL*DSHW^b_SHL
    COND<-STE/(STE+SHE)
    if(SHL>SLM&&COND>=0.95*MTA){
      root<-SST(t)-TTS
    } else if(SHL>SLM&&SST(t)>TTS) {
      root<-COND-0.95*MTA
    } else if(COND>=0.95*MTA&&SST(t)>TTS){
      root<-SHL-SLM
    } else {
      root=1
    }
    yroot<-0-root
    return(yroot)
      
  })
}

spawn_eventfunc <-function(t, y, parms,...){
  with(as.list(c(y,parms)),{
    if(y[6]>365){
      c(y[1],y[2]-y[4]*(1-PSTL),y[3],y[4]*(1-PSTL),0,y[5])
    } else {
      c(y[1],y[2],y[3],y[4],y[5],y[6])
    }
  })
}

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

# dummy data (in lieu of environmental forcings) ######## -------------------------

# construct dummy input data

### temperature  # celcius
SST_mean <- 10
SST_magn <- 10

POC_mean <- 500
POC_magn <- 300

CHL_mean <- 6
CHL_magn  <-5



run_length <- time <- 1:(3*365)



default_input<- data.frame(
  time   = time,

  SST    = SST_mean + (SST_magn*sin(2*pi*(time+180)/365)),
  CHL    = CHL_mean + (CHL_magn*sin(2*pi*(time+180)/365)),
  POC    = POC_mean + (POC_magn*sin(2*pi*(time+180)/365)),
  
  F_in   = 100,
  t_z    = 10,
  D_ext   = 0
)


###################
# parameter load-in: to be moved to run_SF.R eventually

library(rjson)
allparams<-fromJSON(file='shellfish_model_parameters.json')

runparams<-c(unlist(allparams$species$defaults$M_edulis$parameters),harvest_method=0,POC_data=TRUE,ESELORG=23.5)

######################################################
##   run model 
######################################################

run_SF_indiv_model<-function(input,parameters,STE0,SHE0,output='df'){
  #function can be called from R or python, sets up boundary forcing functions from input data and executes the model function, returns a neat data frame of output values

  run_length<-max(input$time)
  times=seq(1,run_length,by=1)
  
  #create boundary forcings
  input_functions = boundary_forcings(input)
  parms = c(parameters, input_functions)

  y0=c(SHE=SHE0,STE=STE0,DSHW=SHE0/(parms$EST*1000),DSTW=STE0/(parms$ECS*1000),spawnday=365,sd2=365)
  print(parms)
  print(y0)
  

    #no harvesting but spawn events
    Out<-ode(times = times,
             func=SF_indiv_model,
             y=y0,
             parms=parms,
             event=list(func=spawn_eventfunc,root=TRUE),
             rootfun=spawn_timed_rootfunc)
              
 
}


######################################################
##   THE MODEL EQUATIONS ###########
######################################################

#Mode based on shellsim equations (Hawkins et al 2013 doi: https://doi.org/10.2983/035.032.0201). Note that the Shelsim model can run on chlorophyll alone or preferably on chrophyll and one of either POC OR POM. Based on the availability of satellite derived POC we implement the model to run either chlrophyll only or chlrophyll and POC. POM input is not implemented.

SF_indiv_model <- function(t,y,parms,...) {

  
  # parameters expected: ####
    # harvest_method
    # Temperature response charateristics
      # CR_max
      # CR_grad
      # CR_inflec
    # m_NSO
    # c_NSO
    # a_NRO
    # b_NRO
    # ON_min
    # ON_max
    # MNEA
    # WCS
    # WCT
    # SCW
    # a_SHL
    # b_SHL
    # PSTW
    # TTS
    # EST
    
  
  
  
  
  # Farm dimensions: 
  
  # Growth parameters: 
   
  
  #time varying inputs: 
   
    #CHL
    #SST
    #POC
    #F_IN
    #T_Z
  
  with(as.list(c(y, parms)), {
    
    #dimensions ####
   
    #V_EFF<-y_farm*(x_farm)*(t_z(t))
    # therefore nutrient change in the outflow of the farm (to depth z+t_z) is scaled to V_MA/V_EFF of the local nutrient change (See dNH4 and dNO3 terms)
    
    #lambda   <- F_in(t)/x_farm
    
    
    if(POC_data){
      SELORG<-CHL(t)*12/(0.38*1000)  #preferentially consumed organic matter in mg/l
      REMORG<-(2.33*POC(t)/1000)-SELORG #Remaining organic matter in mg/l: convert POC in ug/l to POM in mg/l then subtract  SELORG
      EREM<-20.48 # J/mg
    }
    else {
      SELORG<-CHL(t)*50/(0.38*1000) #chlorophyll-rich food in mg/l
      REMORG<-0
      EREM<-0
    }
    

    
    CR<-function(SST.,CR_max,CR_grad,CR_inflec){
      #Clearance rate at temp T, with 3 species specific parameters controlling teh relationship. CR in l/hr/g
      CR_max-CR_grad*(SST.-CR_inflec^2)
    }
    TEF <- CR(SST.=SST(t),CR_max,CR_grad,CR_inflec)/CR(SST.=15,CR_max,CR_grad,CR_inflec) # Temperature effect on feeding
    #NIRSELORG<-ifelse(CHL(t)<0.01,0,((m_NSO*SELORG)+c_NSO)*TEF*DSTW^0.62) # mg/h
    NIRSELORG<-ifelse(CHL(t)<0.01,0,((m_NSO*SELORG))*TEF*(DSTW)^0.62) # mg/h
    #NIRSELORG<-ifelse(CHL(t)<0.01,0,SELORG*CR(SST.=SST(t),CR_max,CR_grad,CR_inflec)*(DSTW)) # mg/h
    #NIRREMORG<-REMORG*DSTW*CR(SST.=SST(t),CR_max,CR_grad,CR_inflec) #mg/h
    
    
    NIRREMORG<-a_NRO*(1-exp(b_NRO*REMORG))*TEF*DSTW^0.62 #mg/h
    
    NEA<-0.8*24*(NIRSELORG*ESELORG + NIRREMORG*0.15*EREM) #J/day
    
    MHL<-4.005*24*(exp(SST(t)*a_TEM)/exp(15*a_TEM))*DSTW^(0.72) # j/day
    THL<-MHL+0.23*NEA #j/day
    #ON<-ON_min+NEA*(ON_max-ON_min)/(MNEA)
    ON<-ON_min+NEA*(ON_max-ON_min)/(MNEA)
    EL<-14*THL/(0.45*ON)
    NEB<-NEA-THL-(EL*0.0248)
    
    COND<-STE/(STE+SHE)
    
    SHG<-ifelse(COND>=MTA,(1-MTA)*NEB,0)#shell growth MTA is mean tissue allocation (between shell and soft tissue, a species specific parameter)
    STG<-ifelse(COND<MTA,NEB,NEB*MTA)#soft tissue growth 
    FW<-SCW*(DSHW*(1+WCS)+DSTW*(1+WCT)) #fresh weight
    DWW<-DSHW*(1+WCS)+DSTW*(1+WCT)
    SHL<-a_SHL*DSHW^b_SHL
    
    
    

    dSHE<-SHG   #shell energy
    dSTE<-STG   #soft tissue energy
    dDSHW<-SHG/(ECS*1000)       #dry shell weight
    dDSTW<-(STG)/(EST*1000) #dry soft tissue weight (g)
    dspawnday<-1
    dsd2<-1
    
    
    
    
    
    output<-list(c(dSHE,dSTE,dDSHW,dDSTW,dspawnday,dsd2),
                 SHG=SHG,
                 STG=STG,
                 COND=COND,
                 SHL=SHL,
                 TEF=TEF,
                 NEB=NEB,
                 EL=EL,
                 NEA=NEA,
                 NIRSELORG=NIRSELORG,
                 NIRREMORG=NIRREMORG,
                 REMORG=REMORG,
                 EREM=EREM,
                 ESELORG=ESELORG,
                 ON=ON,
                 FW=FW,
                 DWW=DWW
                 )
  }) 
}

