# farm post processing to give required paramaters for EMFAF (single value per grid square)
# - kg DW and FW per m2
# - kcal per m2
# - kg protein per m2
# - CO2 balance per m2

# seacarb library allows calculation of carbonate system parameters e.g. alkalinity change due to calcification, effect on air-sea flux
library(seacarb)
library(rjson)



getRunParms<-function(parms_json_file){
  #process json file for run to establish type (shellfish / algae), species and anything else we need 
  parmsJSON<-fromJSON(file=parms_json_file)
  
  #names of macroalgae species modelled
  macroalgae<-c('saccharina','ulva','alaria', 'S_latissima', 'U_lactuta', 'A_esculenta')
  shellfish<-c('M_edulis','C_gigas','P_maximus' )
  
  species<-names(parmsJSON$parameters$species[1])
  if(species %in% macroalgae){
    type<-'macroalgae'
  } else if (species %in% shellfish){
    type<-'shellfish'
  } else {stop('unknown species passed to get_run_parms by json parameter input file')}
  
  species_parms<-get(species,parmsJSON$parameters$species)$parameters
  
  return(c(type=type,species=species,parmsJSON$parameters$farm$default$parameters,species_parms))
}



farmPostProcess_MA<-function(data,json_file,final_only=TRUE){
  #take state variables from MA model output and calculate derived variables as required. Expects a data frame of monthly values
  parms<-getRunParms(json_file)
  with(as.list(c(parms,data)),{
    #calculate DW biomass (B in model)
    DW<-N_f/Q_min #gDW m-3
    DW_line<-DW*h_MA*w_MA / 1000 # kg/m (i.e. per m of rope)
    DW_PUA<-DW*h_MA*density_MA / 1000 # kg/m^2 (multiply be density to account for unused space within farm)
    
  
    
    
    #FW biomass
    FW<-DW/DF_MA #gFW m-3
    FW_line<-DW_line/DF_MA # kg/m (i.e. per m of rope)
    FW_PUA<-DW_PUA/DF_MA # kg/m^2
    
    #Energy
    kcal_PUA<-DW*h_MA*density_MA*kcal_MA # kcal/m^2
    
    #protein
    protein_PUA<-DW_PUA*prot_MA #kg/m^2
    
    #CO2 uptake
    Biomass_CO2<-(N_f/14)*CN_MA*44/1000 #g (CO2) /m^3    (44 is rmm of CO2)
    CO2_uptake_PUA<-Biomass_CO2*h_MA*density_MA/1000 #kg (CO2) / m^2
    
    
    #data_out<-data.frame(DW,DW_line,DW_PUA,FW,FW_line,FW_PUA,kcal_PUA,protein_PUA,Biomass_CO2,CO2_uptake_PUA)
    data_out = list('DW' = DW,
                    'DW_line' = DW_line,
                    'DW_PUA' = DW_PUA,
                    'FW' = FW,
                    'FW_line' = FW_line,
                    'FW_PUA' = FW_PUA,
                    'kcal_PUA' = kcal_PUA,
                    'protein_PUA' = protein_PUA,
                    'Biomass_CO2' = Biomass_CO2,
                    'CO2_uptake_PUA' = CO2_uptake_PUA)
    units_out = list('DW' = 'gDW m-3',
                     'DW_line' = 'kg/m',
                     'DW_PUA' = 'kg/m^2',
                     'FW' = 'gFW m-3',
                     'FW_line' = 'kg/m',
                     'FW_PUA' = 'kg/m^2',
                     'kcal_PUA' = 'kcal/m^2',
                     'protein_PUA' = 'kg/m^2',
                     'Biomass_CO2' = 'g (CO2) /m^3',
                     'CO2_uptake_PUA' = 'kg (CO2) / m^2')
    
    if(final_only){
      for (var in names(data_out)) {
        data_out[[var]] = data_out[[var]][,,dim(data_out[[var]])[3]]
      }
    }
    return(list('data'=data_out, 'units'=units_out))

  })

}

#test_parms_json_file<-'test_macroalgae_model_parameters_input.json'
#test_data<-read.csv('../../../macroalgae/bantry_data/bantry_run.csv',header=TRUE)

#parmsJSON<-fromJSON(file=test_parms_json_file)

