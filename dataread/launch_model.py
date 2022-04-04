import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import STAP
from rpy2.robjects.conversion import localconverter
import pandas as pd
import numpy as np
import datetime
import json
import time
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from types import SimpleNamespace
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from read_netcdf import *


class MA_model:

    def __init__(self, model_path, parms_path):

        # Opens model_path to read the run_MA_model function
        with open(model_path, 'r') as f:
            body = f.read()
            self._model = STAP(body, "run_MA_model")

        ### This section will need an update to use the input json rather than defaults json
        with open(parms_path, 'r') as f:
            parms = json.loads(f.read())

        param_lists = []
        for _,section in parms.items():
            # just take the first option for each section for now
            first_default = section['defaults'][next(iter(section['defaults']))]

            param_lists.append(robjects.ListVector(first_default['parameters']))
            param_lists.append(robjects.ListVector(first_default['options']))

        # get the R concatenator function
        self._conc = robjects.r('c')

        # concatenate all parameter lists
        self._parameters = self._conc(*param_lists)


    def apply_on(self, input_data, latitude, y0=robjects.r('c("NH4" = 0, \
                                                              "NO3" = 0, \
                                                              "N_s" = 1, \
                                                              "N_f" = 1, \
                                                              "D" = 0, \
                                                              "Yield" = 0,\
                                                              "Yield_per_m" = 0)')):
        # Applies the model on the input data (pd.DataFrame), latitude has to be specified

        parms_with_lat = self._conc(self._parameters, robjects.ListVector({"latitude": latitude}))

        # Using a converter between R data.frame and python pd.DataFrame
        with localconverter(robjects.default_converter + pandas2ri.converter):
            out = self._model.run_MA_model(input = input_data,
                                           parameters = parms_with_lat,
                                           y0 = y0)

        return out


class MA_model_scipy:
    def __init__(self, model_path, parms_path):

        # Opens model_path to read the run_MA_model function
        with open(model_path, 'r') as f:
            body = f.read()
            self._model = STAP(body, "MA_model_const")
            self._jacobian = STAP(body, "MA_model_jacobian")

        ### This section will need an update to use the input json rather than defaults json
        with open(parms_path, 'r') as f:
            parms = json.loads(f.read())

        param_lists = []
        self._parameters2 = {}
        for _,section in parms.items():
            # just take the first option for each section for now
            first_default = section['defaults'][next(iter(section['defaults']))]

            param_lists.append(robjects.ListVector(first_default['parameters']))
            param_lists.append(robjects.ListVector(first_default['options']))

            self._parameters2.update(first_default['parameters'])
            self._parameters2.update(first_default['options'])

        # get the R concatenator function
        self._conc = robjects.r('c')
        # concatenate all parameter lists
        self._parameters = self._conc(*param_lists)

        self.names = ["NH4", "NO3", "N_s", "N_f", "D"]

        self.time_reading = 0
        self.time_computing = 0
        self.n_calls = 0
        self.time_reading_J = 0
        self.time_computing_J = 0
        self.n_calls_J = 0


    @staticmethod
    def derivative(t: float, y: np.array, data: pd.DataFrame, latitude: float,
                   startDate: datetime.datetime, self):

        t1 = time.time()
        #yToR = pd.Series(y, index=self.names) # 0.07
        yToR = dict(zip(self.names, y)) # 0.001

        iData = iNearest(t, data['time'], "closest") # 0.04

        
        dataToR = data.iloc[iData] # 0.08
        self.time_reading += time.time() - t1

        t1 = time.time()
        """
        with localconverter(robjects.default_converter + pandas2ri.converter):
            out = self._model.MA_model_const(t = t,
                                             y = yToR,
                                             parms = self._parameters,
                                             data = dataToR,
                                             latitude = latitude
                                             )
        """
        out = self.hadley(t, yToR, dataToR, latitude)
        self.time_computing += time.time() - t1

        self.n_calls += 1

        return out


    def hadley(self, t: float, y, data, latitude: float):
        """Python implementation of the Hadley model adapted by M. Johnson

            y and data must be subscriptable by names (dict, pd.Series,...)
        """

        p = SimpleNamespace(**self._parameters2)

        V_MA = p.y_farm * p.x_farm * p.density_MA * p.h_MA
        V_INT = p.y_farm * p.x_farm * data['t_z']
        V_EFF = p.y_farm * (p.x_farm + data['F_in']) * data['t_z']

        turnover = data['F_in'] / p.x_farm #lambda is taken in python

        Q = p.Q_min * (1 + y['N_s']/y['N_f']) if y['N_f'] > 0 else 0
        B = y['N_f'] / p.Q_min
        g_Q = min(1, (Q - p.Q_min)/(Q - p.K_c))

        T_x = p.T_max if (data['SST'] > p.T_O) else p.T_min
        g_T = np.exp( -2.3 * ((data['SST'] - p.T_O)/(T_x - p.T_O))**2 )

        # Only light scheme == 3 for now
        # still assumes t=0 is jan 1st
        declin = 23.45 * np.cos(np.radians((360/365) * (t + 10)))
        theta = 90 - (latitude + declin)
        sine = np.sin(np.radians(theta))
        I_top = data['PAR'] * np.exp( -data['K_d'] * (p.z - p.h_MA) / sine )
        absorption = (data['K_d']*p.h_MA + y['N_f']*p.a_cs) / sine
        I_av = ( I_top / absorption ) * (1 - np.exp( - absorption ))
        g_E = I_av / ( p.I_s + I_av)

        mu_g_EQT = p.mu * g_E * g_Q * g_T

        f_NH4 = ( p.V_NH4*y['NH4'] / (p.K_NH4 + y['NH4']) ) * ( (p.Q_max - Q) / (p.Q_max - p.Q_min))
        f_NO3 = ( p.V_NO3*y['NO3'] / (p.K_NO3 + y['NO3']) ) * ( (p.Q_max - Q) / (p.Q_max - p.Q_min))

        NH4_removed = f_NH4 * B - p.r_L * y['D'] + p.r_N * y['NH4'] - p.d_m * y['N_s']
        NO3_removed = f_NO3 * B - p.r_N * y['NH4']

        PO4_tot = data['PO4_ext'] * V_EFF / V_MA
        N_to_P_mg = p.N_to_P * 14

        # derivatives
        dNH4 = turnover * (data['NH4_ext'] - y['NH4']) - ( f_NH4*B - p.r_L*y['D'] + p.r_N*y['NH4'] - p.d_m*y['N_s'] ) * V_MA/V_INT

        dNO3 = turnover * (data['NO3_ext'] - y['NO3']) - ( f_NO3*B - p.r_N*y['NH4'] ) * V_MA/V_INT

        dN_s = (f_NH4 + f_NO3)*B - min(mu_g_EQT*y['N_f'], PO4_tot*N_to_P_mg) - p.d_m*y['N_s']

        dN_f = min(mu_g_EQT*y['N_f'], PO4_tot*N_to_P_mg) - p.d_m*y['N_f']

        dD = turnover * (data['D_ext'] - y['D']) + ( p.d_m*y['N_f'] - p.r_L*y['D'] ) * V_MA/V_INT

        return np.array([dNH4, dNO3, dN_s, dN_f, dD])



    @staticmethod
    def jacobian(t: float, y: np.array, data: pd.DataFrame, latitude: float, 
                 startDate: datetime.datetime, self):
        t1 = time.time()
        yToR = dict(zip(self.names, y))

        iData = iNearest(t, data['time'], "closest")
        dataToR = data.iloc[iData]
        self.time_reading_J += time.time() - t1

        t1 = time.time()
        """
        with localconverter(robjects.default_converter + pandas2ri.converter):
            out = self._jacobian.MA_model_jacobian(t = t,
                                                   y = yToR,
                                                   parms = self._parameters,
                                                   data = dataToR,
                                                   latitude = latitude
                                                   )
        """
        out = self.hadley_J(t, yToR, dataToR, latitude)
        self.time_computing_J += time.time() - t1

        self.n_calls_J += 1

        return out

    
    def hadley_J(self, t: float, y, data, latitude: float):
        """Python implementation of the Jacobian of the Hadley model adapted by
        M. Johnson

            y and data must be subscriptable by names (dict, pd.Series,...)
        """

        p = SimpleNamespace(**self._parameters2)

        V_MA = p.y_farm * p.x_farm * p.density_MA * p.h_MA
        V_INT = p.y_farm * p.x_farm * data['t_z']
        V_EFF = p.y_farm * (p.x_farm + data['F_in']) * data['t_z']

        turnover = data['F_in'] / p.x_farm #lambda is taken in python

        Q = p.Q_min * (1 + y['N_s']/y['N_f']) if y['N_f'] > 0 else 0
        B = y['N_f'] / p.Q_min
        g_Q = min(1, (Q - p.Q_min)/(Q - p.K_c))

        T_x = p.T_max if (data['SST'] > p.T_O) else p.T_min
        g_T = np.exp( -2.3 * ((data['SST'] - p.T_O)/(T_x - p.T_O))**2 )

        # Only light scheme == 3 for now
        # still assumes t=0 is jan 1st
        declin = 23.45 * np.cos(np.radians((360/365) * (t + 10)))
        theta = 90 - (latitude + declin)
        sine = np.sin(np.radians(theta))
        I_top = data['PAR'] * np.exp( -data['K_d'] * (p.z - p.h_MA) / sine )
        absorption = (data['K_d']*p.h_MA + y['N_f']*p.a_cs) / sine
        I_av = ( I_top / absorption ) * (1 - np.exp( - absorption ))
        g_E = I_av / ( p.I_s + I_av)

        f_NH4 = ( p.V_NH4*y['NH4'] / (p.K_NH4 + y['NH4']) ) * ( (p.Q_max - Q) / (p.Q_max - p.Q_min))
        f_NO3 = ( p.V_NO3*y['NO3'] / (p.K_NO3 + y['NO3']) ) * ( (p.Q_max - Q) / (p.Q_max - p.Q_min))

        PO4_tot = data['PO4_ext'] * V_EFF / V_MA
        N_to_P_mg = p.N_to_P * 14

        ### Useful derivatives and values
        df_NH4_dNH4 = ( (p.V_NH4*p.K_NH4) / (p.K_NH4+y['NH4'])**2 ) * ((p.Q_max-Q)/(p.Q_max-p.Q_min))
        df_NO3_dNO3 = ( (p.V_NO3*p.K_NO3) / (p.K_NO3+y['NO3'])**2 ) * ((p.Q_max-Q)/(p.Q_max-p.Q_min))

        df_NH4_dQ = ((p.V_NH4*y['NH4'])/(p.K_NH4+y['NH4'])) * (-1/(p.Q_max-p.Q_min))
        df_NO3_dQ = ((p.V_NO3*y['NO3'])/(p.K_NO3+y['NO3'])) * (-1/(p.Q_max-p.Q_min))

        dQ_dNs = p.Q_min / y['N_f']
        dQ_dNf = - p.Q_min * y['N_s'] / y['N_f']**2

        dB_dNf = 1 / p.Q_min

        dg_dQ = g_T * g_E * (p.Q_min - p.K_c)/(Q - p.K_c)**2

        g = g_T * g_E * g_Q

        limited_P = (p.mu * g * y['N_f']) > (PO4_tot * N_to_P_mg)

        # All Jacobian values
        dJ_NH4_dNH4 = - turnover + ( - B * df_NH4_dNH4 - p.r_N ) * V_MA/V_INT
        dJ_NH4_dNO3 = 0
        dJ_NH4_dNs = (p.d_m - B * df_NH4_dQ * dQ_dNs ) * V_MA/V_INT
        dJ_NH4_dNf = ( - dB_dNf * f_NH4 - B * df_NH4_dQ * dQ_dNf ) * V_MA/V_INT
        dJ_NH4_dD = p.r_L * V_MA/V_INT

        dJ_NO3_dNH4 = p.r_N * V_MA/V_INT
        dJ_NO3_dNO3 = - turnover + ( - B * df_NO3_dNO3 )  * V_MA/V_INT
        dJ_NO3_dNs = ( - B * df_NO3_dQ * dQ_dNs ) * V_MA/V_INT
        dJ_NO3_dNf = ( - dB_dNf * f_NO3 - B * df_NO3_dQ * dQ_dNf ) * V_MA/V_INT
        dJ_NO3_dD = 0

        dJ_Ns_dNH4 = B * df_NH4_dNH4
        dJ_Ns_dNO3 = B * df_NO3_dNO3
        dJ_Ns_dNs = ( df_NH4_dQ + df_NO3_dQ ) * dQ_dNs * B - (0 if limited_P else p.mu * dg_dQ * dQ_dNs * y['N_f']) - p.d_m
        dJ_Ns_dNf = ( df_NH4_dQ + df_NO3_dQ ) * dQ_dNf * B - (0 if limited_P else p.mu * dg_dQ * dQ_dNf * y['N_f'] + p.mu * g )
        dJ_Ns_dD = 0

        dJ_Nf_dNH4 = 0
        dJ_Nf_dNO3 = 0
        dJ_Nf_dNs = (0 if limited_P else p.mu * dg_dQ * dQ_dNs * y['N_f'])
        dJ_Nf_dNf = (0 if limited_P else p.mu * dg_dQ * dQ_dNf * y['N_f'] + p.mu * g ) - p.d_m
        dJ_Nf_dD = 0

        dJ_D_dNH4 = 0
        dJ_D_dNO3 = 0
        dJ_D_dNs = 0
        dJ_D_dNf = p.d_m * V_MA/V_INT
        dJ_D_dD = - turnover - p.r_L * V_MA/V_INT

        output = np.array([[dJ_NH4_dNH4, dJ_NO3_dNH4, dJ_Ns_dNH4, dJ_Nf_dNH4, dJ_D_dNH4],
                           [dJ_NH4_dNO3, dJ_NO3_dNO3, dJ_Ns_dNO3, dJ_Nf_dNO3, dJ_D_dNO3],
                           [dJ_NH4_dNs,  dJ_NO3_dNs,  dJ_Ns_dNs,  dJ_Nf_dNs,  dJ_D_dNs],
                           [dJ_NH4_dNf,  dJ_NO3_dNf,  dJ_Ns_dNf,  dJ_Nf_dNf,  dJ_D_dNf],
                           [dJ_NH4_dD,   dJ_NO3_dD,   dJ_Ns_dD,   dJ_Nf_dD,   dJ_D_dD]
                            ]).T
        return output



if __name__ =="__main__":

    # Read from CSV for now
    bantryData = pd.read_csv('bantry_data/bantry_3m.csv', delimiter=';')
    bantryData['date'] = [datetime.datetime.fromisoformat(str_date) for str_date in bantryData['date']]

    dataToR = pd.DataFrame({
        'time': [(date - bantryData['date'][0]).days + 1 for date in bantryData['date']],
        'SST': 15,
        'PAR': bantryData['par'],
        'NH4_ext': bantryData['Ammonium'],
        'NO3_ext': bantryData['Nitrate'],
        'PO4_ext': 50,
        'K_d': 0.1,
        'F_in': 100,
        'h_z_SML': 30,
        't_z': 10,
        'D_ext': 0.1
    })

    dataToR = pd.DataFrame({
        'time': [0, 1],
        'SST': 15,
        'PAR': bantryData['par'][100],
        'NH4_ext': bantryData['Ammonium'][100],
        'NO3_ext': bantryData['Nitrate'][100],
        'PO4_ext': 50,
        'K_d': 0.1,
        'F_in': 100,
        'h_z_SML': 30,
        't_z': 10,
        'D_ext': 0.1
    })

    y0 = pd.DataFrame({"NH4": [0],
                       "NO3": 0,
                       "N_s": 1,
                       "N_f": 1,
                       "D": 0,
                       "Yield": 0,
                       "Yield_per_m": 0})

    model = MA_model("macroalgae_model.R", "macroalgae_model_parameters.json")

    print(dataToR)

    out = model.apply_on(dataToR, 51., y0)

    print(out)


    zone = "IBI"

    mainpath = '/media/share/data_merged/'

    #dataRef = pd.read_csv('/profils/qjutard/shellfish_and_algae-MODEL/dataimport/src/dataCmd.csv', delimiter=';')
    dataRef = pd.read_csv('./dataCmd.csv', delimiter=';')

    ### Initialize the netcdf reading interface

    #paramNames = ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current', 'ocean_mixed_layer_thickness', 'par']
    paramNames = ['Ammonium', 'Nitrate', 'Temperature', 'northward_Water_current', 'eastward_Water_current']
    fileNames = [mainpath + f"{zone}/{param}/{zone}_{param}_merged.nc" for param in paramNames]

    #dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone)][0] for param in paramNames]
    dataRows = [dataRef.index[(dataRef['Parameter']==param) & (dataRef['Place']==zone)][-1] for param in paramNames]

    variableNames = dataRef.iloc[dataRows]['variable'].tolist()
    latitudeNames = dataRef.iloc[dataRows]['latName'].fillna('latitude').tolist()
    longitudeNames = dataRef.iloc[dataRows]['longName'].fillna('longitude').tolist()
    timeNames = dataRef.iloc[dataRows]['timeName'].fillna('time').tolist()
    depthNames = dataRef.iloc[dataRows]['depthName'].fillna('depth').tolist()
    unitConversions = dataRef.iloc[dataRows]['unitFactor'].fillna(1).tolist()

    algaeData = AllData(fileNameList=fileNames,
                        parameterNameList=paramNames,
                        variableNameList=variableNames,
                        latitudeNameList=latitudeNames,
                        longitudeNameList=longitudeNames,
                        timeNameList=timeNames,
                        depthNameList=depthNames,
                        unitConversionList=unitConversions
    )

    startDate = datetime.datetime(2021, 1, 1, 12)
    endDate = datetime.datetime(2022, 1, 1, 12)

    model = MA_model_scipy("macroalgae_model.R", "macroalgae_model_parameters.json")

    input_data = algaeData.getTimeSeries(51.5874, -9.8971, (startDate, endDate), 3)
    data = pd.DataFrame({
                'time': [(date - input_data['date'][0]).days for date in input_data['date']],
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
    print(data)

    y0 = np.array([0, 0, 0, 1000, 0])

    #dy0 = MA_model_scipy.derivative(y0, 0, data, 49, startDate, model)
    #J0 = MA_model_scipy.jacobian(y0, 0, data, 49, startDate, model)
    #print(dy0)
    #print(J0)

    model = MA_model_scipy("macroalgae_model.R", "macroalgae_model_parameters.json")

    n_runs = 1
    t0 = time.time()

    for _ in range(n_runs):
        result = solve_ivp(MA_model_scipy.derivative, (0, 365), y0, args=(data, 51.5874, startDate, model),
                           rtol=0.05, atol=[0.7, 0.7, 500, 500, 0.05], method='BDF')

    print(f"One run: {(time.time() - t0)/n_runs}\n" +
          f"Reading: {model.time_reading}\n" +
          f"Computing: {model.time_computing}\n" +
          f"N of calls: {model.n_calls} ({(model.time_reading+model.time_computing)/(model.n_calls or 1)} per call)\n"
          )
    
    #print(result[:, 3])


    t0 = time.time()

    model = MA_model_scipy("macroalgae_model.R", "macroalgae_model_parameters.json")

    result_J = solve_ivp(MA_model_scipy.derivative, (0, 365), y0, args=(data, 51.5874, startDate, model),
                         jac=MA_model_scipy.jacobian,
                         rtol=0.05, atol=[0.7, 0.7, 100, 100, 0.5], method='BDF')

    print(f"One run with jacobian: {time.time() - t0}\n" +
          f"Reading: {model.time_reading}\n" +
          f"Computing: {model.time_computing}\n" +
          f"N of calls: {model.n_calls} ({(model.time_reading+model.time_computing)/model.n_calls} per call)\n" +
          f"Reading jacobian: {model.time_reading_J}\n" +
          f"Computing jacobian: {model.time_computing_J}\n" +
          f"N of calls jacobian: {model.n_calls_J} ({(model.time_reading_J+model.time_computing_J)/(model.n_calls_J or 1)} per call)"
           )

    plt.plot(result.t, result.y[3,:])
    plt.plot(result_J.t, result_J.y[3,:])
    plt.savefig('/media/share/results/foo.png')

    # Used to detect possible errors in the Jacobian
    """
    h = 1
    y_last = result_J.y[:,-1]
    J_last = MA_model_scipy.jacobian(365, y_last, data, 51.5874, startDate, model)
    f_last = MA_model_scipy.derivative(365, y_last, data, 51.5874, startDate, model)
    #print(np.linalg.eig(J_last))
    print(J_last)

    errors = np.zeros((5,5))
    for i in range(5):
        for j in range(5):
            e_j = np.zeros((5))
            e_j[j] = 1

            f_last_h = MA_model_scipy.derivative(365, y_last+h*e_j, data, 51.5874, startDate, model)

            errors[i,j] = ( (1/h)*(f_last_h - f_last) - J_last[:,j] )[i]
    print(errors)
    print(abs(errors) > abs(J_last))
    """