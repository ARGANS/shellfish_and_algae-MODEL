# farm post processing to give required paramaters for EMFAF (single value per grid square)
# - kg DW and FW per m2
# - kcal per m2
# - kg protein per m2
# - CO2 balance per m2

# seacarb library allows calculation of carbonate system parameters e.g. alkalinity change due to calcification, effect on air-sea flux
library(seacarb)
library(rjson)

#names of macroalgae species modelled
macroalgae<-c('saccharina','ulva','alaria', 'S_latissima', 'U_lactuta', 'A_esculenta')
shellfish<-c('M_edulis','C_gigas','P_maximus' )

getRunParms<-function(parms_json_file){
  #process json file for run to establish type (shellfish / algae), species and anything else we need 
  parmsJSON<-fromJSON(file=parms_json_file)
  species<-names(parmsJSON$parameters$species[1])
  if(species %in% macroalgae){
    type<-'macroalgae'
  } else if (species %in% shellfish){
    type<-'shellfish'
  } else stop('unknown species passed to get_run_parms by json parameter input file')
  return(c(type=type,species=species,parmsJSON$parameters$farm$default$parameters))
}

farmPostProcess_protein<-function(data){
 #return the protein content yield in kg per m2 across the whole farm area
  data$protein_PUA[nrow(data)]
}

farmPostProcess_kcal<-function(data){
  #return the energy yield in kcal per m2 across the whole farm area
  data$kcal_PUA[nrow(data)]
}

farmPostProcess_CO2<-function(data){
  #simple estimate of farm CO2 uptake in tonnes CO2 per farm - need refining to account for buffering by carbonate system
  tonnesC<-cumsum(data$Farm_CO2_demand)[nrow(data)]/1e6
  tonnesC*3.67
}

farmPostProcess_DW_yield<-function(data,json_file){
  #DW yield in kg/m2
  parms<-getRunParms(json_file)
  if(parms$type=='macroalgae'){
    data$B_vol*parms$h_ma
  }
    else if(parms$type=='shellfish'){
      
  } else stop ('could not parse farm type in farmPostProcess function')
}

test_parms_json_file<-'../../../macroalgae/macroalgae_model_parameters_input.json'
test_data<-read.csv('../../../macroalgae/bantry_data/bantry_run.csv',header=TRUE)

parmsJSON<-fromJSON(file=test_parms_json_file)

