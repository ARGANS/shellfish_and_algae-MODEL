#!/usr/bin/env Rscript

library(ncdf4)

source('farmPostprocess/farmprocess.R')

# Get the command line parameters
args = commandArgs(trailingOnly=TRUE)

input_file = args[1]
parms_file = args[2]

nc = nc_open(input_file, write=T)

state_vars = c('NH4', 'NO3', 'N_s', 'N_f', 'D')

data = list()
py_default_fillval = 9969209968386869046778552952102584320
for (varName in state_vars) {
    #data[[varName]] = ncvar_get(nc, varName, start=c(1, 1, 12), count=c(-1, -1, 1))
    data[[varName]] = ncvar_get(nc, varName)
    data[[varName]][data[[varName]]==py_default_fillval] = NaN
}

res = farmPostProcess_MA(data, parms_file, TRUE)


for (interest_var in names(res$data)) {
    new_var = ncvar_def(name = interest_var,
                        units = res$units[[interest_var]],
                        #dim = list(nc$dim$latitude, nc$dim$longitude), 
                        dim = list(nc$dim$longitude, nc$dim$latitude), 
                        missval = py_default_fillval)
    nc = ncvar_add(nc, new_var)
    ncvar_put(nc, new_var, res$data[[interest_var]])
}

nc_close(nc)