import numpy as np
import matplotlib.pyplot as plt
from toymodel_update import toymodel_update
from water_energy_variables import water_energy_variables
from weighted_average_outputs import weighted_average_outputs
import aggregate as agg

# 1 - Load Data
Q = np.genfromtxt('Data/Q_mm.csv')  # Q in mm
Dates = Q[:, 0:3]
Precip = np.genfromtxt('Data/Precip.csv')  # P in mm
e = np.genfromtxt('Data/e.csv')  # e in kPa
u2 = np.genfromtxt('Data/u2.csv')  # Wind Speed at 2m in m/s
S_in = np.genfromtxt('Data/S_in.csv')  # Solar radiation in W/m2
Temp = np.genfromtxt('Data/Temp.csv')  # Temp in Celsius
Lat_Lon_A_Z = np.genfromtxt('Data/Lat_Lon_Area_Z.csv')  # Elevation

Q = Q[:, 3:]
Precip = Precip[:, 3:]
e = e[:, 3:]
u2 = u2[:, 3:]
S_in = S_in[:, 3:]
Temp = Temp[:, 3:]

# 2.1 - Generate Daily PET: Calculate Net Longwave
Ra, Rso = Rso_calc(Dates, Lat_Lon_A_Z[:, 1], Lat_Lon_A_Z[:, 2], Lat_Lon_A_Z[:, 4])  # Assuming Rso_calc is defined
Rso = Rso / (1e-6 * 60 * 60 * 24)
Ra = Ra / (1e-6 * 60 * 60 * 24)

f = 1.35 * (S_in / Rso) - 0.35
e_prime = 0.34 - 0.14 * np.sqrt(e / 1000)
L_net = -f * e_prime * 5.67e-8 * (Temp + 273.15)**4

# 2.2 - Generate Daily PET
cp = 1.1013e3
lambda_v = 1e3 * (2.501 - 0.002361 * Temp)
rho_air = 1.23
ra = 208 / u2
rs = 70
Rn_RC = S_in * (1 - 0.23) + L_net
e_sat = 0.6108 * np.exp((17.27 * Temp) / (237.3 + Temp))
D = e_sat - e
Delta = 4098 * e_sat / (237.3 + Temp)**2
Z = Lat_Lon_A_Z[:, 4]
Press = 101.3 * ((293 - 0.0065 * Z) / 293)**5.26
gamma = cp * Press / (0.622 * lambda_v * 1e3)

Wm2_mm = (lambda_v * 1e3)**-1 * 86400

E_RC = (Delta * Rn_RC + rho_air * cp * D / ra) / (Delta + gamma * (1 + rs / ra))
E_RC_mm = E_RC * Wm2_mm
Rn_RC_mm = Rn_RC * Wm2_mm

# 3.1 Mean annual analysis
n_years = len(np.unique(Dates[:, 0]))
PET_MEAN_MONTHLY, PET_MEAN_ANNUAL = agg.make_means_new(E_RC_mm, Dates)  # Assuming make_means_new is defined
Q_MEAN_MONTHLY, Q_MEAN_ANNUAL = agg.make_means_new(Q, Dates)
P_MEAN_MONTHLY, P_MEAN_ANNUAL = agg.make_means_new(Precip, Dates)
Temp_MEAN_MONTHLY, Temp_MEAN_ANNUAL = agg.make_means(Temp, Dates)  # Assuming make_means is defined

PHI = PET_MEAN_ANNUAL / P_MEAN_ANNUAL
E_P = (P_MEAN_ANNUAL - Q_MEAN_ANNUAL) / P_MEAN_ANNUAL

# Budyko plot
AI_plot_lin = np.arange(0, 50, 0.01)
Budyko = lambda AI: (AI * (1 - np.exp(-AI)) * np.tanh(1 / AI))**0.5

# 3.2 Run Hydrologic Model (assuming toymodel_update and other necessary functions are defined)
pick_catchment = 10  # Example of picking a catchment

# Selecting input data for the specific catchment
INPUT = np.column_stack((Precip[:, pick_catchment], 
                         E_RC_mm[:, pick_catchment],
                         e[:, pick_catchment],
                         u2[:, pick_catchment],
                         S_in[:, pick_catchment],
                         Temp[:, pick_catchment],
                         L_net[:, pick_catchment]))

Dates_sim = Dates.copy()

# Setting up variables for different environments
vars_forest = water_energy_variables(1)  # Assuming this function is defined
vars_grass = water_energy_variables(2)

# Parameter setup
PAR = [50, 150, 30, 1]  # Example parameter values

# Setting specific parameters and conditions
vars_1 = vars_forest.copy()
vars_1['Patm'] = Press[pick_catchment]
vars_1['g0'] = 10
vars_1['z0'] = 0.826
vars_1['a'] = 0.12

Su_max = PAR[1]
vars_1['FC'] = 0.35 * Su_max
vars_1['WP'] = 0.11 * Su_max

# Running the hydrological model simulation
OUT_1 = toymodel_update(INPUT, PAR, vars_1)
OUT_1 = agg.aggregate(OUT_1, Dates_sim)  # Assuming aggregate function is defined

# Calculating additional metrics
PHI_SIM = OUT_1['PET_mean_annual'] / OUT_1['P_mean_annual']
EF_SIM = (OUT_1['P_mean_annual'] - OUT_1['QT_mean_annual']) / OUT_1['P_mean_annual']

# Assuming OUT_1 and other variables are dictionaries or similar structures
Q_SIM_mean_monthly = OUT_1['QT_mean_mon']
Q_SIM = OUT_1['QT']
plot_from = 1000
plot_to = 1300

# Simulation 2
vars_2A = vars_forest.copy()
vars_2A['Patm'] = Press[pick_catchment]
vars_2A['g0'] = 10
vars_2A['z0'] = 0.83
vars_2A['a'] = 0.12
vars_2A['FC'] = 0.38 * Su_max
vars_2A['WP'] = 0.11 * Su_max

vars_2B = vars_grass.copy()
vars_2B['Patm'] = Press[pick_catchment]
vars_2B['g0'] = 3
vars_2B['z0'] = 0.24
vars_2B['a'] = 0.25
vars_2B['FC'] = 0.38 * Su_max
vars_2B['WP'] = 0.11 * Su_max

PAR2B = PAR.copy()
PAR2B[0] = 25

weight_A = 0.3
weight_B = 0.7

OUT_2A = toymodel_update(INPUT, PAR, vars_2A)
OUT_2B = toymodel_update(INPUT, PAR2B, vars_2B)

OUT_2 = weighted_average_outputs(OUT_2A, OUT_2B, weight_A, weight_B)
OUT_2 = agg.aggregate(OUT_2, Dates_sim)

# Produce final figures
plt.figure(1)
plt.plot(Ra[plot_from:plot_to, pick_catchment], '.b')
plt.plot(Rso[plot_from:plot_to, pick_catchment], '.k')
plt.plot(S_in[plot_from:plot_to, pick_catchment], '.-r')
plt.legend(['Extraterrestrial solar radiation', 'Clear sky solar radiation', 'S in'])
plt.ylabel('W/m^2')
plt.xlabel('days')

plt.figure(2)
plt.subplot(1, 2, 1)
plt.plot(PHI, E_P, 'ok', linewidth=1)
plt.plot(PHI[pick_catchment], E_P[pick_catchment], '.r', markersize=20)
plt.plot(AI_plot_lin, Budyko(AI_plot_lin), 'k-', linewidth=1)
plt.legend(['Catchments', 'Your pick', 'Budyko'])
plt.ylabel('E/P')
plt.xlabel('Phi')
plt.axis([0, 5, 0, 1])

plt.subplot(1, 2, 2)
plt.plot(range(1, 13), P_MEAN_MONTHLY[:, pick_catchment], '-b', linewidth=1)
plt.plot(range(1, 13), PET_MEAN_MONTHLY[:, pick_catchment], '-r', linewidth=1)
plt.ylabel('mm')
plt.xlabel('Month')
plt.twinx()
plt.plot(range(1, 13), Temp_MEAN_MONTHLY[:, pick_catchment], '--m', linewidth=1)
plt.ylabel('deg. C')
plt.legend(['P', 'PET', 'Temp'])

# Plot_gSM is a function to plot soil moisture, assumed to be defined
Plot_gSM(vars_1['FC'], vars_1['WP'], Su_max)

plt.figure(4)
plt.subplot(1, 2, 1)
plt.plot(PHI, E_P, 'ok', linewidth=1)
plt.plot(PHI[pick_catchment], E_P[pick_catchment], '.r', markersize=20)
plt.plot(PHI_SIM, EF_SIM, '.b', markersize=20)
plt.plot(AI_plot_lin, Budyko(AI_plot_lin), 'k-', linewidth=1)
plt.legend(['Catchments', 'Your catchment', 'Your catchment (simulated)', 'Budyko'])
plt.ylabel('E/P')
plt.xlabel('Phi')
plt.axis([0, 5, 0, 1])

plt.subplot(1, 2, 2)
plt.plot(range(1, 13), Q_SIM_mean_monthly, '--k', linewidth=1)
plt.plot(range(1, 13), Q_MEAN_MONTHLY[:, pick_catchment], '--r', linewidth=1)
plt.legend(['Simulated', 'Observed'])

plt.figure(5)
plt.subplot(1, 2, 1)
plt.plot(Q_SIM[plot_from:plot_to], '-k', linewidth=0.5)
plt.plot(Q[plot_from:plot_to, pick_catchment], '-r', linewidth=0.5)
plt.legend(['Simulated', 'Observed'])

plt.subplot(1, 2, 2)
plt.plot(OUT_1['Su'], '-k', linewidth=0.5)
plt.legend(['Simulated'])
plt.twinx()
plt.plot(OUT_1['QS'], '-r', linewidth=0.5)

