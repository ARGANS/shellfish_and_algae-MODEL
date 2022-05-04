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
    DSHW=SHE/(ECS*1000)
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
    if(y[4]>365){
      c(y[1],y[2]*(1-PSTL),0,y[3],y[5],y[6],y[7],y[8])
    } else {
      c(y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8])
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




###################
# parameter load-in: to be moved to run_SF.R eventually

library(rjson)
allparams<-fromJSON(file='shellfish_model_parameters.json')

runparams<-c(unlist(allparams$species$defaults$M_edulis$parameters),harvest_method=0,POC_data=TRUE,PHYC_data=FALSE,ESELORG=23.5,Bioavailable_detrital_fraction=0.15,Assimilation_efficiency=0.8,y_farm=1000,x_farm=1000,density_SF_farm=0.3)

######################################################
##   run model 
######################################################

run_SF_model<-function(input,parameters,y0,output='df'){
  #function can be called from R or python, sets up boundary forcing functions from input data and executes the model function, returns a neat data frame of output values

  run_length<-max(input$time)
  times=seq(1,run_length,by=1)
  
  #create boundary forcings
  input_functions = boundary_forcings(input)
  parms = c(parameters, input_functions)


  
  if(parms['harvest_method']==0){
    #no harvesting
    Out<-ode(times = times,
             func=SF_model,
             y=y0,
             parms=parms,
             event=list(func=spawn_eventfunc,root=TRUE),
             rootfun=spawn_timed_rootfunc,
             maxsteps=5e4,atol=1e-5)
              
  } else if(parms['harvest_method']==1){
    
    #do some harvesting
    
    
    Out_timed_harvest <- ode(times = times,
                             func = SA_model,
                             y = y0,
                             parms = parms,
                             event=list(func=harvest_eventfunc, root=TRUE),
                             rootfun=harvest_timed_rootfunc)
    Out<-Out_timed_harvest
    
  } 
  
  if(output=='df'){
    Out.<-as.data.frame(Out)
    colnames(Out.)[1]<-c('model_time')
    cbind(input,Out.)
  } else {Out}
}


######################################################
##   THE MODEL EQUATIONS ###########
######################################################

#Mode based on shellsim equations (Hawkins et al 2013 doi: https://doi.org/10.2983/035.032.0201). Note that the Shelsim model can run on chlorophyll alone or preserable on chrophyll and one of either POC OR POM. Based on the availability of satellite derived POC we implement the model to run either chlrophyll only or chlrophyll and POC. POM input is not implemented.

SF_model <- function(t,y,parms,...) {
 
  
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
   
    #CHL or PHYC (ug/l)
    #SST
    #POC
    #F_IN
    #T_Z
  
  with(as.list(c(y, parms)), {
    
    #dimensions ####
    #farm dimensions are defined by input parms x_farm, y_farm (both 1000m i.e. grid square dimensions by default) and z_farm (depth to bottom of macroalgae)
    #Shellfish are grown on notional 'lines' running the full width of y_farm (perpendicular to notional flow in x direction). The 'spacing' of the lines is determined by input parm density_MA. The total volume occupied by macroalgae is calculated as y_farm*x_farm*density_SF*h_SF, where h_SF is a run-specific parameter representing the height of the aquaculture system e.g trestles vs long-line etc
    V_SF<-y_farm*x_farm*density_SF_farm*h_SF
    # The volume of water that the V_SF is interacting with due to vertical mixing is 
    V_INT<-y_farm*x_farm*t_z(t)
    # therefore concentration changes in the outflow of the farm (to depth t_z) is scaled to V_MA/V_INT of the local nutrient change (See dNH4, dCHL, dPOC terms)
    
    # within V_SF we need to calculate the number of individuals so we can scale the outputs of the individual-based model appropriately
    #n_SF<-V_SF*stocking_density #stocking density in individuals per m3 (a species_specific parmeter)
    
    
    lambda   <- (F_in(t)/x_farm) #
    
    # convenience terms for calculating outputs
    
    per_farm <- POP*V_SF
    per_unit_area <-per_farm/(x_farm*y_farm)
    
 ## the equations below represent an individual shellfish, representative of the average shellfish accross the whole farm. The farm-scale state variable equations at the end of the model and the derived output parameters scale to the farm level using the above terms.    
    if(POC_data){
      if(PHYC_data){
      SELORG<-PHYC_farm*2/1000
      } else {
      SELORG<-CHL_farm*12/(0.38*1000)  #preferentially consumed organic matter in mg/l
      }
      REMORG<-max(0,(2.33*POC_farm/1000)-SELORG) #Remaining organic matter in mg/l: convert POC in ug/l to POM in mg/l then subtract  SELORG
      EREM<-20.48 # J/mg
    }
    else {
      SELORG<-CHL_farm*50/(0.38*1000) #chlorophyll-rich food in mg/l
      REMORG<-0
      EREM<-0
    }
    
    DSHW<-SHE/(ECS*1000)       #dry shell weight (g)
    DSTW<-(STE)/(EST*1000) #dry soft tissue weight (g)
    
    CR<-function(SST.,CR_max,CR_grad,CR_inflec){
      #Clearance rate at temp T, with 3 species specific parameters controlling teh relationship. CR in l/hr/g
      CR_max-CR_grad*(SST.-CR_inflec^2)
    }
   TEF <- CR(SST.=SST(t),CR_max,CR_grad,CR_inflec)/CR(SST.=15,CR_max,CR_grad,CR_inflec) # Temperature effect on feeding
    #NIRSELORG<-ifelse(CHL(t)<0.01,0,((m_NSO*SELORG)+c_NSO)*TEF*DSTW^0.62) # mg/h
    NIRSELORG<-ifelse(CHL(t)<0.01,0,((m_NSO*SELORG))*TEF*(DSTW)^0.62) # mg/h
    #NIRSELORG<-ifelse(CHL(t)<0.01,0,SELORG*CR(SST.=SST(t),CR_max,CR_grad,CR_inflec)*(DSTW^0.62)) # mg/h
    #NIRREMORG<-REMORG*(DSTW^0.62)*CR(SST.=SST(t),CR_max,CR_grad,CR_inflec) #mg/h

    
    NIRREMORG<-a_NRO*(1-exp(b_NRO*REMORG))*TEF*DSTW^0.62 #mg/h
    
    NEA<-Assimilation_efficiency*24*(NIRSELORG*ESELORG + NIRREMORG*Bioavailable_detrital_fraction*EREM) #J/day
    
    MHL<-4.005*24*(exp(SST(t)*a_TEM)/exp(15*a_TEM))*DSTW^(0.72) # j/day
    THL<-MHL+0.23*NEA #j/day
    #ON<-ON_min+NEA*(ON_max-ON_min)/(MNEA)
    ON<-ON_min+NEA*(ON_max-ON_min)/(MNEA)
    EL<-14*THL/(0.45*ON) # ug NH4-N excreted (note this N is produced arbitrarily by the model - not linked to N content of input but to estimated Oxygen consumption and O:N ratio given conditions of organism...)
    NEB<-NEA-THL-(EL*0.0248)
    
    
    COND<-STE/(STE+SHE)
    
    SHG<-ifelse(COND>=MTA&&NEB>0,(1-MTA)*NEB,0)#shell growth MTA is mean tissue allocation (between shell and soft tissue, a species specific parameter)
    if(NEB>0){
      STG<-ifelse(COND<MTA,NEB,NEB*MTA)#soft tissue growth
      starvation_mortality=0
      }else{
        STG=0
        starvation_mortality=POP*NEB/STE # calculate individuals lost as population multiplied by the ratio of individual NEB to soft tissue energy
      }   
    FW<-SCW*(DSHW*(1+WCS)+DSTW*(1+WCT)) #fresh weight
    DWW<-DSHW*(1+WCS)+DSTW*(1+WCT)
    SHL<-a_SHL*DSHW^b_SHL
    
    
    #calculate quantity of POC and CHL taken up by individual
    
    
      
    if(POC_data){ 
      PHYC_uptake<-NIRSELORG*24*1000/2 #ug per day per individual
      CHL_uptake<-NIRSELORG*(0.38*1000)*24/12 #ug per day per individual
      POC_uptake<-CHL_uptake+(NIRREMORG*Bioavailable_detrital_fraction)*1000*24/2.33 #ug per day per individual
      } else {
      CHL_uptake<-NIRSELORG*(0.38*1000)*24/50 #ug per day per individual
      POC_uptake<-0 
    }
    
#individual shellfish state variables
    
    dSHE<-SHG   #shell energy
    dSTE<-STG   #soft tissue energy
    #dDSHW<-SHG/(ECS*1000)       #dry shell weight (g)
    #dDSTW<-(STG)/(EST*1000) #dry soft tissue weight (g)
    

#spawning logic state variables 
    dspawnday<-1
    dsd2<-1
    
#farm level state variables
    
    if(POC_data){
      dPOC_farm<-lambda*(POC(t)-POC_farm)-(POC_uptake*POP*V_SF/V_INT)/1000
    } else {
      dPOC_farm<-0
    }
    if(PHYC_data){
      dCHL_farm<-0
      dPHYC<-lambda*(PHYC(t)-PHYC_farm)-(PHYC_uptake*POP*V_SF/V_INT)/1000 #ug/l/day
    }else{
      dCHL_farm<-lambda*(CHL(t)-CHL_farm)-(CHL_uptake*POP*V_SF/V_INT)/1000 #divide by 1000 to convert to ug/L/day
      dPHYC_farm<-0
    }
    #population numbers per m3
   dPOP<-starvation_mortality#-POP*(specific_mortality+0.00001)#for some reason specific mortality values of 0.01,0.001,0.0001... cause integrator to fail - weird numerical error - add a tiny bit to fix it)
    
    output<-list(c(dSHE,dSTE,dspawnday,dsd2,dCHL_farm,dPOC_farm,dPOP,dPHYC_farm),
                 #individual data
                 DSTW=DSTW,
                 DSHW=DSHW,
                 SHG=SHG,
                 STG=STG,
                 COND=COND,
                 SHL=SHL,
                 #TEF=TEF,
                 NEB=NEB,
                 EL=EL,
                 NEA=NEA,
                 SELORG=SELORG,
                 NIRSELORG=NIRSELORG,
                 NIRREMORG=NIRREMORG,
                 REMORG=REMORG,
                 EREM=EREM,
                 ESELORG=ESELORG,
                 ON=ON,
                 FW=FW,
                 DWW=DWW,
                 CHL_uptake=CHL_uptake,
                 POC_uptake=POC_uptake,
                 #volume terms
                 V_SF=V_SF,
                 V_INT=V_INT,
                 #output terms
                 NH4_production=EL*per_farm/1e6,#g NH4-N per day  per farm
                 DWW_total=DWW*per_farm/1e6, #tonnes per farm
                 DWW_PUA=DWW*per_unit_area/1e3, #kg per m2
                 FW_total=FW*per_farm/1e6,#tonnes per farm
                 FW_PUA=FW*per_unit_area/1e3, #kg per m2
                 DW_total=(DSTW+DSHW)*per_farm/1e6,#tonnes per farm
                 DW_PUA=(DSTW+DSHW)*per_unit_area/1e3, #kg per m2
                 DSTW_total=(DSTW)*per_farm/1e6,#tonnes per farm
                 DSTW_PUA=(DSTW)*per_unit_area/1e3, #kg per m2
                 POC_uptake_total=POC_uptake*per_farm/1e6,#g per day  per farm
                 POC_uptake_PUA=POC_uptake*per_unit_area/1e6 #g per m2 per day
                 
                 )
  }) 
}

#X<-run_SF_model(default_input,runparams,STE0 = 1000, SHE0 = 100)
