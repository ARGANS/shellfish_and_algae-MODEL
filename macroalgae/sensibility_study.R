
require(reshape2)
require(RColorBrewer)

source('run_MA.R')

sensibility_study <- function(param_name, values, default_parms, default_input_data) {
  # param: string with the name of the parameter that will vary
  # values: vector of values that param should take
  # default_parms: named list of parameters for the model
  # input_data: data_frame of inputs to the model, it also defines the run duration
  
  if (param_name %in% names(default_parms)) {
    input_data = default_input_data
    boundary_forcings(input_data)
  } else if (param_name %in% names(default_input_data)) {
    parms = default_parms
  } else {
    print(paste0(param_name, " is not a parameter in the model/inputs."))
    return(1)
  }
  
  y0 = c(NH4=NH4_in(1), NO3=NO3_in(1), N_s=0.01, N_f=0.01, D=0.1)
  
  param_dependency = data.frame(NULL)
  for (param_val in values) {
    
    if (param_name %in% names(default_parms)) {
      parms = default_parms
      parms[param_name] = param_val
    } else {
      input_data = default_input_data
      input_data[param_name] = param_val
      boundary_forcings(input_data)
    }
    
    simulated = ode(times=input_data$time, func=MA_model, y=y0, parms=parms)
    
    param_simulated = data.frame(param_val, simulated)
    names(param_simulated)[1] = param_name
    
    param_dependency = rbind(param_dependency, param_simulated)
  }
  
  return(param_dependency)
}


ab = sensibility_study("F_in", c(0.25, 1), default_parms, dummy_input)
#ab = sensibility_study("mu", c(0, 0.33), default_parms, dummy_input)

plot_sens_time = ggplot(data=ab, mapping=aes(x=time, y=NH4, color=as.factor(F_in))) +
  geom_point() +
  scale_color_discrete() +
  theme_bw()

plot_sens_time






