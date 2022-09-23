import pandas as pd
import numpy as np
import datetime
import json
import time
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from types import SimpleNamespace
#import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class MA_model_scipy:
    def __init__(self, parms: dict):

        self._parameters = {}
        for _,section in parms.items():
            # take the first option for each section which should be the only one
            first_default = next(iter(section.values()))

            self._parameters.update(first_default['parameters'])
            self._parameters.update(first_default['options'])

        self.names = ["NH4", "NO3", "N_s", "N_f", "D"]

        self.time_reading = 0
        self.time_computing = 0
        self.n_calls = 0
        self.time_reading_J = 0
        self.time_computing_J = 0
        self.n_calls_J = 0


    @staticmethod
    def derivative(t: float, y: np.array, data:dict, latitude: float, self):
        """
        data is a dict of functions giving the value of the data at time t.
        """

        yToR = dict(zip(self.names, y)) # 0.001
        #t1 = time.time()
        dataToR = {name:fun(t) for name,fun in data.items()}
        #self.time_reading += time.time() - t1
        #t1 = time.time()
        out = self.hadley(t, yToR, dataToR, latitude)
        #self.time_computing += time.time() - t1

        #self.n_calls += 1

        return out
    
    @staticmethod
    def derivative_fast(t: float, y: np.array, data:dict, latitude: float, self):
        """
        data is a dict of functions giving the value of the data at time t.

        same as derivative(). Used when data can be passed directly to hadley()
        """

        yToR = dict(zip(self.names, y)) # 0.001

        out = self.hadley(t, yToR, data, latitude)

        return out


    def hadley(self, t: float, y, data, latitude: float):
        """Python implementation of the Hadley model adapted by M. Johnson

            y and data must be subscriptable by names (dict, pd.Series,...)
        """

        p = SimpleNamespace(**self._parameters)

        V_MA = p.y_farm * p.x_farm * p.density_MA * p.h_MA
        V_INT = p.y_farm * p.x_farm * data['t_z']

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
        if theta > 0 :
            sine = np.sin(np.radians(theta))
            I_top = data['PAR'] * np.exp( -data['K_d'] * (p.z - p.h_MA) / sine )
            absorption = (data['K_d']*p.h_MA + y['N_f']*p.a_cs) / sine
            I_av = ( I_top / absorption ) * (1 - np.exp( - absorption ))
            g_E = I_av / ( p.I_s + I_av)
        else:
            g_E = 0

        mu_g_EQT = p.mu * g_E * g_Q * g_T

        f_NH4 = ( p.V_NH4*y['NH4'] / (p.K_NH4 + y['NH4']) ) * ( (p.Q_max - Q) / (p.Q_max - p.Q_min)) if y['NH4'] > 0 else 0
        f_NO3 = ( p.V_NO3*y['NO3'] / (p.K_NO3 + y['NO3']) ) * ( (p.Q_max - Q) / (p.Q_max - p.Q_min)) if y['NO3'] > 0 else 0

        #NH4_removed = f_NH4 * B - p.r_L * y['D'] + p.r_N * y['NH4'] - p.d_m * y['N_s']
        #NO3_removed = f_NO3 * B - p.r_N * y['NH4']

        PO4_tot = data['PO4_ext'] * V_INT / V_MA
        N_to_P_mg = p.N_to_P * 14

        # derivatives
        dNH4 = turnover * (data['NH4_ext'] - y['NH4']) - ( f_NH4*B - p.r_L*y['D'] + p.r_N*y['NH4'] - p.d_m*y['N_s'] ) * V_MA/V_INT

        dNO3 = turnover * (data['NO3_ext'] - y['NO3']) - ( f_NO3*B - p.r_N*y['NH4'] ) * V_MA/V_INT

        dN_s = (f_NH4 + f_NO3)*B - min(mu_g_EQT*y['N_f'], PO4_tot*N_to_P_mg) - p.d_m*y['N_s']

        dN_f = min(mu_g_EQT*y['N_f'], PO4_tot*N_to_P_mg) - p.d_m*y['N_f']

        dD = turnover * (data['D_ext'] - y['D']) + ( p.d_m*y['N_f'] - p.r_L*y['D'] ) * V_MA/V_INT

        return np.array([dNH4, dNO3, dN_s, dN_f, dD])

    @staticmethod
    def derivative_fast_advection(t: float, y: np.array, data: dict,cNO3, cNH4, advNO3, advNH4, nbrx, nbry, dt, latitude: float, self):
        """
        data is a dict of functions giving the value of the data at time t.

        same as derivative(). Used when data can be passed directly to hadley()
        """

        # yToR = dict(zip(self.names, y)) # 0.001

        out = self.hadley_advection(t, y, data, cNO3, cNH4, advNO3, advNH4, nbrx, nbry, dt, latitude)

        return out

    def hadley_advection(self, t: float, y, data, dt, latitude: float, nitrification_NH4=None):
        """Python implementation of the Hadley model adapted by M. Johnson

            y and data must be subscriptable by names (dict, pd.Series,...)
        """

        p = SimpleNamespace(**self._parameters)

        V_MA = p.y_farm * p.x_farm * p.density_MA * p.h_MA
        V_INT = p.y_farm * p.x_farm * data['t_z']

        #turnover = data['F_in'] / p.x_farm  # lambda is taken in python
        Q = (p.Q_min * (1 + y['N_s'] / np.maximum(y['N_f'],1e-6))) * (y['N_f'] > 0)
        B = y['N_f'] / p.Q_min
        g_Q = np.minimum(1, (Q - p.Q_min) / (Q - p.K_c))

        T_x = p.T_max * (data['SST'] > p.T_O) + p.T_min * (1 - (data['SST'] > p.T_O))
        g_T = np.exp(-2.3 * ((data['SST'] - p.T_O) / (T_x - p.T_O)) ** 2)

        # Only light scheme == 3 for now
        # still assumes t=0 is jan 1st
        declin = 23.45 * np.cos(np.radians((360 / 365) * (t + 10)))
        theta = 90 - (latitude + declin)
        sine = np.sin(np.radians(theta))*(theta > 0)+1*(theta <= 0)
        I_top = data['PAR'] * np.exp(-data['K_d'] * (p.z - p.h_MA) / sine)
        absorption = (data['K_d'] * p.h_MA + y['N_f'] * p.a_cs) / sine
        I_av = (I_top / absorption) * (1 - np.exp(- absorption))
        g_E = (I_av / (p.I_s + I_av))*(theta > 0)

        mu_g_EQT = p.mu * g_E * g_Q * g_T

        f_NH4 = (p.V_NH4 * y['NH4'] / (p.K_NH4 + y['NH4'])) * ((p.Q_max - Q) / (p.Q_max - p.Q_min))
        f_NO3 = (p.V_NO3 * y['NO3'] / (p.K_NO3 + y['NO3'])) * ((p.Q_max - Q) / (p.Q_max - p.Q_min))
        f_NH4[y['NH4'] < 0] = 0
        f_NO3[y['NO3'] < 0] = 0

        #NH4_removed = f_NH4 * B - p.r_L * y['D'] + p.r_N * y['NH4'] - p.d_m * y['N_s']
        #NO3_removed = f_NO3 * B - p.r_N * y['NH4']

        PO4_tot = data['PO4_ext'] * V_INT / V_MA
        N_to_P_mg = p.N_to_P * 14

        if nitrification_NH4 is None:
            nitrification_term = p.r_N * y['NH4']
        else:
            nitrification_term = p.r_N * nitrification_NH4

        # derivatives
        dNH4 = ( - f_NH4 * B + p.r_L * y['D'] - nitrification_term + p.d_m * y['N_s']) * V_MA / V_INT

        dNO3 = (- f_NO3 * B + nitrification_term) * V_MA / V_INT

        dN_s = (f_NH4 + f_NO3) * B - np.minimum(mu_g_EQT * y['N_f'], PO4_tot * N_to_P_mg) - p.d_m * y['N_s']

        dN_f = np.minimum(mu_g_EQT * y['N_f'], PO4_tot * N_to_P_mg) - p.d_m * y['N_f']

        dD = (p.d_m * y['N_f'] - p.r_L * y['D']) * V_MA / V_INT

        return np.array([dNH4, dNO3, dN_s, dN_f, dD])

    @staticmethod
    def jacobian(t: float, y: np.array, data: dict, latitude: float, self):
        #t1 = time.time()
        yToR = dict(zip(self.names, y))

        dataToR = {name:fun(t) for name,fun in data.items()}
        #self.time_reading_J += time.time() - t1

        #t1 = time.time()
        out = self.hadley_J(t, yToR, dataToR, latitude)
        #self.time_computing_J += time.time() - t1

        #self.n_calls_J += 1

        return out

    @staticmethod
    def jacobian_fast(t: float, y: np.array, data: dict, latitude: float, self):
        """Same as jacobian when data can be directly passed tp hadley_J
        """
        yToR = dict(zip(self.names, y))

        out = self.hadley_J(t, yToR, data, latitude)

        return out

    
    def hadley_J(self, t: float, y, data, latitude: float):
        """Python implementation of the Jacobian of the Hadley model adapted by
        M. Johnson

            y and data must be subscriptable by names (dict, pd.Series,...)
        """

        p = SimpleNamespace(**self._parameters)

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
        if theta > 0 :
            sine = np.sin(np.radians(theta))
            I_top = data['PAR'] * np.exp( -data['K_d'] * (p.z - p.h_MA) / sine )
            absorption = (data['K_d']*p.h_MA + y['N_f']*p.a_cs) / sine
            I_av = ( I_top / absorption ) * (1 - np.exp( - absorption ))
            g_E = I_av / ( p.I_s + I_av)
        else:
            g_E = 0

        f_NH4 = ( p.V_NH4*y['NH4'] / (p.K_NH4 + y['NH4']) ) * ( (p.Q_max - Q) / (p.Q_max - p.Q_min))
        f_NO3 = ( p.V_NO3*y['NO3'] / (p.K_NO3 + y['NO3']) ) * ( (p.Q_max - Q) / (p.Q_max - p.Q_min))

        PO4_tot = data['PO4_ext'] * V_EFF / V_MA
        N_to_P_mg = p.N_to_P * 14

        ### Useful derivatives and values
        df_NH4_dNH4 = ( (p.V_NH4*p.K_NH4) / (p.K_NH4+y['NH4'])**2 ) * ((p.Q_max-Q)/(p.Q_max-p.Q_min)) if y['NH4'] > 0 else 0
        df_NO3_dNO3 = ( (p.V_NO3*p.K_NO3) / (p.K_NO3+y['NO3'])**2 ) * ((p.Q_max-Q)/(p.Q_max-p.Q_min)) if y['NO3'] > 0 else 0

        df_NH4_dQ = ((p.V_NH4*y['NH4'])/(p.K_NH4+y['NH4'])) * (-1/(p.Q_max-p.Q_min)) if y['NH4'] > 0 else 0
        df_NO3_dQ = ((p.V_NO3*y['NO3'])/(p.K_NO3+y['NO3'])) * (-1/(p.Q_max-p.Q_min)) if y['NO3'] > 0 else 0

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
    pass