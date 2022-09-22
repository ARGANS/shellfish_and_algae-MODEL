from pickle import TRUE
import pandas as pd
import numpy as np
import datetime
import json
import time
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from types import SimpleNamespace
import sys, os
#import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class SF_model_scipy:
    def __init__(self, parms: dict):
        self._parameters = {}
        for _,section in parms.items():
            # take the first option for each section which should be the only one
            first_default = next(iter(section.values()))
            self._parameters.update(first_default["parameters"])
            self._parameters.update(first_default["options"])

        self.names = ["SHE", "STE", "spawnday", "sd2", "POP","CHL","FW","DWW","SHL","NH4_production","CO2_production"]
        self.time_reading = 0
        self.time_computing = 0
        self.n_calls = 0
        self.time_reading_J = 0
        self.time_computing_J = 0
        self.n_calls_J = 0

    @staticmethod
    def SF_derivative(y, data, self):
        p = SimpleNamespace(**self._parameters)
        y = spawn_tined_root(y, data)
        V_SF = p.y_farm*p.x_farm*p.density_SF_farm*p.h_SF
        V_INT = p.y_farm*p.x_farm*p.h_SF*1.4 #Zmix
        turnover = data['F_in'] / p.x_farm
        per_farm = y["POP"]*V_SF
        per_unit_area = per_farm/(p.x_farm*p.y_farm)

        if p.POC_data == TRUE :
            if p.PHYC_data == TRUE : 
                SELORG = data['PHYC']*2/1000
            else:  
                SELORG = data['CHL_ext']*12/(0.38*1000)
                REMORG = max(0, (2.33*data['POC']/1000) - SELORG)
                EREM = 20.48
        else : 
            SELORG = data['CHL_ext']*50/0.38
            REMORG = 0
            EREM = 0
        DSHW = y['SHE']/(p.ECS*1000)
        DSTW = y['STE']/(p.EST*1000)
        TEF = CR(data['SST'], p.CR_max,p.CR_grad, p.CR_inflec)/CR(15, p.CR_max,p.CR_grad, p.CR_inflec)
        NIRSELORG = (p.m_NSO*SELORG)*TEF*(DSTW**0.62)*(data['CHL_ext'] > 0.01)
        NIRREMORG = p.a_NRO*(1- np.exp(p.b_NRO*REMORG))*TEF*(DSTW**0.62) 
        NEA = p.Assimilation_efficiency*24*(NIRSELORG*p.ESELORG + NIRREMORG*p.Bioavailable_detrital_fraction*EREM) 
        MHL = 4.005*24*(np.exp(data['SST']*p.a_TEM)/np.exp(15*p.a_TEM))*DSTW**(0.72)
        THL = MHL+0.23*NEA 
        ON = p.ON_min+NEA*(p.ON_max-p.ON_min)/(p.MNEA)
        EL = 14*THL/(0.45*ON)
        ExcretaCNratio=6
        ELC = EL*(12/14)*ExcretaCNratio
        NEB = NEA-THL-(EL*0.0248)
        COND = y['STE']/(y['STE']+y['SHE'])
        SHG = ((1-p.MTA)*NEB)*(COND >= p.MTA)*(NEB > 0) + 0*(NEB <= 0)*(COND < p.MTA)
        STG = 0*(NEB <= 0) + NEB*(COND < p.MTA)*(NEB > 0) + NEB*p.MTA*(COND >= p.MTA)*(NEB > 0)
        starvation_mortality = y["POP"]*NEB/y['STE']*(NEB < 0) + 0*(NEB >= 0)
        FW = p.SCW*(DSHW*(1+p.WCS)+DSTW*(1+p.WCT)) 
        DWW = DSHW*(1+p.WCS)+DSTW*(1+p.WCT)
        SHL = p.a_SHL*DSHW**p.b_SHL
        if p.POC_data == TRUE : 
            PHYC_uptake = NIRSELORG*24*1000/2 
            CHL_uptake = NIRSELORG*(0.38*1000)*24/12
            POC_uptake = CHL_uptake+(NIRREMORG*p.Bioavailable_detrital_fraction)*1000*24/2.33
        else :
            CHL_uptake = NIRSELORG*(0.38*1000)*24/50
            POC_uptake = 0 
        #individual shellfish state variables
        dSHE = SHG   #shell energy
        dSTE = STG   #soft tissue energy
        #spawning logic state variables 
        dspawnday = 1
        dsd2 = 1
        #farm level state variables
        if p.POC_data == TRUE: 
            dPOC_farm = turnover*(y['POC'] - data['POC'])-(POC_uptake*y["POP"]*V_SF/V_INT)/1000
        else:
            dPOC_farm = 0
        if p.PHYC_data == TRUE:
            dCHL_farm = 0 
            dPHYC_farm = turnover*(y['PHYC'] - data['PHYC'])-(PHYC_uptake*y["POP"]*V_SF/V_INT)/1000 #ug/l/day
        else: 
            dCHL_farm = turnover*(y['CHL'] - data['CHL_ext'])-(CHL_uptake*y["POP"]*V_SF/V_INT)/1000 #divide by 1000 to convert to ug/L/day
            dPHYC_farm = 0

        dPOP = starvation_mortality

        NH4_production=EL*per_farm/1e6,#g NH4-N per day  per farm
        CO2_production=ELC*per_farm/1e6, #g CO2-C per day per farm
        DWW_total=DWW*per_farm/1e6, #tonnes per farm
        DWW_PUA=DWW*per_unit_area/1e3, #kg per m2
        FW_total=FW*per_farm/1e6,#tonnes per farm
        FW_PUA=FW*per_unit_area/1e3, #kg per m2
        DW_total=(DSTW+DSHW)*per_farm/1e6,#tonnes per farm
        DW_PUA=(DSTW+DSHW)*per_unit_area/1e3, #kg per m2
        DSTW_total=(DSTW)*per_farm/1e6,#tonnes per farm
        DSTW_PUA=(DSTW)*per_unit_area/1e3, #kg per m2
        POC_uptake_total=POC_uptake*per_farm/1e6,#g per day  per farm
        POC_uptake_PUA=POC_uptake*per_unit_area/1e6, #g per m2 per day
        STE_total=y["STE"]*per_farm/1e6, #MJ soft tissue energy per farm
        STE_PUA=y["STE"]*per_unit_area/1e6, #MJ soft tissue energy per unit area
        kcal_PUA=y["STE"]*per_unit_area/4184, #kcal sfot tissue per unit area
        protein_PUA=DSTW*per_unit_area*p.protein_DW_fraction/1e3 #kg/m2 '''      

        return  np.array([dSHE, dSTE, dspawnday, dsd2, dPOP, dCHL_farm, dPOC_farm, dPHYC_farm, FW, DWW, SHL, NH4_production, CO2_production, kcal_PUA, protein_PUA])

def CR(sst, cr_max, cr_gard, cr_inflec):
     #Clearance rate at temp T, 
     # with 3 species specific parameters controlling the relationship. CR in l/hr/g
    cr = cr = cr_max-cr_gard*(sst-cr_inflec)**2
    CR = cr*(cr>0) +0*(cr<=0)
    return CR     

def spawn_evenfunc(y,self):
    p = SimpleNamespace(**self._parameters)
    print(y)
    y = dict(zip(["CHL", "SHE", "STE", "spawnday", "sd2", "POP"], [y["CHL"],y["SHE"]*(1-p.PSTL)*(y["sd2"]> 365 )+ (y["SHE"]*(1-p.PSTL)*(y["sd2"] < 365 )),0,y["spawnday"],y["sd2"],y["POP"]])) #reduce dstw soft tissu weight
    print(y)
    return y

def spawn_tined_root(y, data, self):
    p = SimpleNamespace(**self._parameters)
    DSHW = y['SHE']/(p.ECS*1000)
    SHL = p.a_SHL*(DSHW**p.b_SHL)
    COND = y['STE']/(y['STE']+y['SHE'])
    root = (data['SST'] - p.TTS)*(SHL > p.SLM)*(COND >= 0.95*p.MTA) + (COND - 0.95*p.MTA)*(SHL > p.SLM)*(data['SST'] > p.TTS) + (SHL - p.SLM)*(COND >= 0.95*p.MTA)*(data['SST'] > p.TTS) + 1*(SHL < p.SLM)*(COND < 0.95*p.MTA)*(data['SST'] < p.TTS)
    yroot = 0-root
    return yroot

'''
def run_SF_model(input,parameters,y0,output='df'):
      #function can be called from R or python, sets up boundary forcing functions from input data and executes the model function, returns a neat data frame of output values
    run_length<-max(input$time)
    times=seq(1,run_length,by=1)

    #create boundary forcings
    input_functions = boundary_forcings(input)
    parms = c(parameters, input_functions)
    Out<-ode(times = times,
                func=SF_model,
                y=y0,
                parms=parms,
                event=list(func=spawn_eventfunc,root=TRUE),
                rootfun=spawn_timed_rootfunc,
                maxsteps=5e4,atol=1e-5)
                
    if(output=='df'){
    Out.<-as.data.frame(Out)
    colnames(Out.)[1]<-c('model_time')
    cbind(input,Out.)
    } else {Out}
    }
    return 0

'''


if __name__ =="__main__":
    zone = "IBI"
    mainpath = os.path.join( '..', 'Data') #'/media/share/data_merged/'
    dataRef = pd.read_csv('./dataCmd.csv', delimiter=';')

    ### Initialize the netcdf reading interface

    paramNames = ['Temperature', 'Chlorophyll_a', 'northward_Water_current', 'eastward_Water_current']
    fileNames = [mainpath + f"{zone}/{param}/{zone}_{param}_merged.nc" for param in paramNames]

    dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone)][0] for param in paramNames]
    #dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone)][-1] for param in paramNames]

    variableNames = dataRef.iloc[dataRows]['variable'].tolist()
    latitudeNames = dataRef.iloc[dataRows]['latName'].fillna('latitude').tolist()
    longitudeNames = dataRef.iloc[dataRows]['longName'].fillna('longitude').tolist()
    timeNames = dataRef.iloc[dataRows]['timeName'].fillna('time').tolist()
    depthNames = dataRef.iloc[dataRows]['depthName'].fillna('depth').tolist()
    unitConversions = dataRef.iloc[dataRows]['unitFactor'].fillna(1).tolist()

    ShellfishData = AllData(fileNameList=fileNames,
                        parameterNameList=paramNames,
                        variableNameList=variableNames,
                        latitudeNameList=latitudeNames,
                        longitudeNameList=longitudeNames,
                        timeNameList=timeNames,
                        depthNameList=depthNames,
                        unitConversionList=unitConversions
    )

    startDate = datetime.datetime(2020, 1, 1, 12)
    endDate = datetime.datetime(2020, 1, 1, 12)

    input_data = ShellfishData.getTimeSeries(51.5874, -9.8971, (startDate, endDate), 3)

    time_axis = [(date - input_data['date'][0]).days for date in input_data['date']]
    interpKind = "nearest"
    data_fun = {
        'SST': interp1d(time_axis, input_data['Temperature'], kind=interpKind, assume_sorted=True),
        'CHL-a': interp1d(time_axis, input_data['Chlorophyll-a'], kind=interpKind, assume_sorted=True),
        'F_in': interp1d(time_axis, np.sqrt(input_data['northward_Water_current']**2 + input_data['eastward_Water_current']**2), kind=interpKind, assume_sorted=True),
        'h_z_SML': lambda _ : 30,
        't_z': lambda _ : 10,
    }

    y0 = np.array([0, 0, 0, 1000., 0])

    model = SF_model_scipy("shellfish_model_parameters.json")

    n_runs = 1
    t0 = time.time()
    print(t0)

    for _ in range(n_runs):
        result_J = solve_ivp(SF_model_scipy.derivative, (0, 364), y0, args=(data_fun, 51.5874, model),
                            jac=MA_model_scipy.jacobian, t_eval=time_axis,
                            rtol=0.05, method='BDF')
    if not result_J.success:
        print(result_J.message)

    print(f"One run with jacobian: {(time.time() - t0)/n_runs}\n" +
          f"Reading: {(model.time_reading)/n_runs}\n" +
          f"Computing: {(model.time_computing)/n_runs}\n" +
          f"N of calls: {model.n_calls/n_runs} ({(model.time_reading+model.time_computing)/(model.n_calls or 1)} per call)\n" +
          f"Reading jacobian: {model.time_reading_J/n_runs}\n" +
          f"Computing jacobian: {model.time_computing_J/n_runs}\n" +
          f"N of calls jacobian: {model.n_calls_J/n_runs} ({(model.time_reading_J+model.time_computing_J)/(model.n_calls_J or 1)} per call)"
           )

    t_m = np.array([0])
    result_m = np.expand_dims(y0, 1)
    initTime = datetime.datetime(2021, 1, 1, 0)
    t0 = time.time()
    for _ in range(n_runs):
        for month in range(1,13):
            #print(f"MONTH: {month}")
            startTime = datetime.datetime(2021, month, 1, 0)
            if month == 12:
                endTime = datetime.datetime(2022, 1, 1, 0)
            else:
                endTime = datetime.datetime(2021, month+1, 1, 0)

            data, dims = algaeData.getData(longitude = -9.8971,
                                       latitude = 51.5874,
                                       depth = 3,
                                       time = (startTime, endTime),
                                       averagingDims = ('time',),
                                       weighted = False
                                        )
            data_in = {
                    'SST': np.average(data['Temperature']),
                    'PAR': 500,
                    'NH4_ext': np.average(data['Ammonium']),
                    'NO3_ext': np.average(data['Nitrate']),
                    'PO4_ext': 50,
                    'K_d': 0.1,
                    'F_in': np.average(np.sqrt(data['northward_Water_current']**2 + data['eastward_Water_current']**2)),
                    'h_z_SML': 30,
                    't_z': 10,
                    'D_ext': 0.1
                }
            
            days_start = (startTime - initTime).days
            days_end = (endTime - initTime).days
            out_m = solve_ivp(MA_model_scipy.derivative_fast, (days_start, days_end), y0, args=(data_in, 51.5874, model),
                                jac=MA_model_scipy.jacobian_fast, #t_eval=[n_days],
                                rtol=0.05, method='BDF')

            result_m = np.concatenate((result_m, out_m.y), 1)
            t_m = np.concatenate((t_m, out_m.t))

            y0 = np.squeeze(out_m.y[:,-1])

    print(f"One run, monthly averaged: {(time.time() - t0)/n_runs}\n")

    data_csv = pd.DataFrame({
        'time': time_axis,
        'SST': input_data['Temperature'],
        'PAR': 500,
        'NH4_ext': input_data['Ammonium'],
        'NO3_ext': input_data['Nitrate'],
        'PO4_ext': 50,
        'K_d': 0.1,
        'F_in': np.sqrt(input_data['northward_Water_current']**2 + input_data['eastward_Water_current']**2),
        'h_z_SML': 30,
        't_z': 10,
        'D_ext': 0.1
    })

    print(data_csv)
    print(model._parameters)
    print(pd.DataFrame(result_J.y.T, columns=model.names))

    data_csv.to_csv("/media/share/results/example/input_data.csv", sep=';')
    pd.DataFrame(result_J.y.T, columns=model.names).to_csv("/media/share/results/example/output_data.csv", sep=';')
    with open('/media/share/results/example/input_parms.json', 'w') as outfile:
        json.dump(model._parameters, outfile)