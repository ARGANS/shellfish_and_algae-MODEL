require(ggplot2)

source('run_MA.R')
source('hadley_model_test.R')

sensitivity_study <- function(param_name, values, default_parms, default_input_data,
                              model_func, y0) {
  # param: string with the name of the parameter that will vary
  # values: vector of values that param should take
  # default_parms: named list of parameters for the model
  # input_data: data_frame of inputs to the model, it also defines the simulation duration
  # model_func: the function of the model on which the sens study is run
  
  
  if (param_name %in% names(default_parms)) {
    input_data = default_input_data
    boundary_forcings(input_data)
  } else if (param_name %in% names(default_input_data)) {
    parms = default_parms
  } else {
    print(paste0(param_name, " is not a parameter in the model/inputs."))
    return(1)
  }
  
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
    
    simulated = ode(times=input_data$time, func=model_func, y=y0, parms=parms)
    
    param_simulated = data.frame(param_val, simulated)
    names(param_simulated)[1] = param_name
    
    param_dependency = rbind(param_dependency, param_simulated)
  }
  
  return(param_dependency)
}


parms_hadley = c(parms_porphyra, parms_farm)
boundary_forcings(input_hadley)
y0 = c(NH4=NH4_ref(1), NO3=NO3_ref(1), N_s=0.01, N_f=1, D=0.1)

# run the sensibility study
sens_hadley = sensitivity_study("F_in", c(0.05*3, 1*3), parms_hadley, input_hadley, 
                                model=Hadley_model_as_published, y0=y0)

# Compute total N as in figure 5 of the paper
sens_hadley$total_N = (sens_hadley$N_s + sens_hadley$N_f) * parms_hadley['h_MA']

plot_sens_time = ggplot(data=sens_hadley, mapping=aes(x=time, y=total_N, color=as.factor(F_in))) +
  geom_point() +
  scale_color_discrete() +
  theme_bw()

# plot the first year
plot_sens_time %+% subset(sens_hadley, time<=365)




