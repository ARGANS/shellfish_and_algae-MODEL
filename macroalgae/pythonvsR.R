

### test run with bantry bay data - MJ 1-/02/2022

source('run_MA.R')
source('calc_par.R')


library('rjson')
#the new data will already have unit conversion done..


input<-read.csv('python vs R comparison/input_data.csv',sep=';')

output<-read.csv('python vs R comparison/output_data.csv',sep=';')
parms<-fromJSON(file = 'python vs R comparison/input_parms.json')
parms<-c(parms,CN_MA=12,prot_MA=8,kcal_MA=2)


test<-run_MA_model(input=input,
                         parameters=parms,
                         y0=c(NH4=0,NO3=0,N_s=0,N_f=10000,D=0,Yield_farm=0, Yield_per_m=0)
                         )



plot(test$NH4)
points(output$NH4[1:366],pch=0.2,col='red')
plot(test$NO3)
lines(output$NO3[1:366])
plot(test$N_f)
lines(output$N_f[1:366])
plot(test$N_s)
lines(output$N_s[1:366])
plot(test$D)
lines(output$D[1:366])


