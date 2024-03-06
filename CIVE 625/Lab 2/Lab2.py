import numpy as np
import matplotlib.pyplot as plt
from Rso_calc import Rso_calc
import toymodel as tm

# 1 - Load Data
Q = np.loadtxt('CIVE 625/Lab 2/Data/Q_mm.csv', delimiter=',')   # Q in mm
Dates = Q[:, :3]
Precip = np.loadtxt('CIVE 625/Lab 2/Data/Precip.csv', delimiter=',') # P in mm
e = np.loadtxt('CIVE 625/Lab 2/Data/e.csv', delimiter=',')      # e in kPa
u2 = np.loadtxt('CIVE 625/Lab 2/Data/u2.csv', delimiter=',')    # Wind Speed at 2m in m/s
S_in = np.loadtxt('CIVE 625/Lab 2/Data/S_in.csv', delimiter=',') # Solar radiation in W/m2
Temp = np.loadtxt('CIVE 625/Lab 2/Data/Temp.csv', delimiter=',') # Temp in Celsius
Lat_Lon_A_Z = np.loadtxt('CIVE 625/Lab 2/Data/Lat_Lon_Area_Z.csv', delimiter=',') # Lat, Lon, Area, Z

# Remove first 3 columns from each array
Q = Q[:, 3:]
Precip = Precip[:, 3:]
e = e[:, 3:]
u2 = u2[:, 3:]
S_in = S_in[:, 3:]
Temp = Temp[:, 3:]

# 2.1 - Generate Daily PET: Calculate Net Longwave
# Assuming Rso_calc is defined elsewhere
Ra, Rso = Rso_calc(Dates, Lat_Lon_A_Z[:, 0], Lat_Lon_A_Z[:, 1], Lat_Lon_A_Z[:, 3])

Rso = Rso / (10**-6 * 60 * 60 * 24) # Convert to W/m2
Ra = Ra / (10**-6 * 60 * 60 * 24)   # Convert to W/m2

from_index = 100
to_index = 300
chose = 14  # Python indexing starts from 0, so 15 in MATLAB is 14 in Python

plt.figure(1)
plt.clf()
plt.plot(Ra[from_index:to_index, chose], '.b', label='Extraterrestrial solar radiation')
plt.plot(Rso[from_index:to_index, chose], '.k', label='Clear sky solar radiation')
plt.plot(S_in[from_index:to_index, chose], '.-r', label='S in')
plt.legend()
plt.ylabel('W/m^2')
plt.xlabel('days')

f = S_in / Rso
e_prime = 0.34 - 0.14 * np.sqrt(e / 1000)
L_net = -f * e_prime * 5.67 * 10**-8 * (Temp + 273.15)**4

# 2.2 - Generate Daily PET:
cp = 1.1013  # specific heat at constant pressure for air kJ/(kg.K)
lambda_ = 1000 * (2.501 - 0.002361 * Temp)  # Latent Heat of Vaporization (kJ/kg)
rho_air = 1.23
ra = 208 / u2
rs = 70
Rn_RC = S_in * (1 - 0.23) + L_net  # Net Radiation over a Reference Crop, in W/m2
e_sat = 0.6108 * np.exp((17.27 * Temp) / (237.3 + Temp))  # Saturated Vapor Pressure in kPa
D = e_sat - e  # Vapor Pressure Deficit in kPa
Delta = 4098 * e_sat / (237.3 + Temp)**2  # Slope of esat versus temp curve (kPa/C)
Z = Lat_Lon_A_Z[:, 3]  # Elevation (m)
Press = 101.3 * ((293 - 0.0065 * Z) / 293)**5.26  # Pressure as function of elevation (kPa)
gamma = cp * Press / (0.622 * lambda_)  # Psychrometric Constant (kPa/degC)

Wm2_mm = (lambda_ * 1000)**-1 * 86400  # Conversion Factor from W/m2 to mm/day

E_RC = (Delta * Rn_RC + rho_air * cp * D / ra) / (Delta + gamma * (1 + rs / ra))

E_RC_mm = E_RC * Wm2_mm
Rn_RC_mm = Rn_RC * Wm2_mm

# Plot Example
from_index = 500
to_index = 800
chose = 12  # Python indexing starts from 0

plt.figure(1)
plt.clf()
plt.subplot(1, 2, 1)
plt.plot(E_RC[from_index:to_index, chose], '.b', label='Ref. Crop ET')
plt.plot(Rn_RC[from_index:to_index, chose], '.-k', label='Net Radiation')
plt.legend()
plt.ylabel('W/m2')
plt.xlabel('days')

plt.subplot(1, 2, 2)
plt.plot(E_RC_mm[from_index:to_index, chose], '.b', label='Ref. Crop ET')
plt.plot(Rn_RC_mm[from_index:to_index, chose], '.-k', label='Net Radiation')
plt.legend()
plt.ylabel('mm')
plt.xlabel('days')

# Run Hydrologic Model
pick_catchment = 7  # Python indexing starts from 0

INPUT = np.zeros((Precip.shape[0], 2))
INPUT[:, 0] = Precip[:, pick_catchment]
INPUT[:, 1] = E_RC_mm[:, pick_catchment]

PAR = [50, 450, 5, 1]  # Parameters

# Assuming toymodel is defined elsewhere
Ea, QF, R, QS, QT, Sf, Su, Ss, St, AL, IE, SE = tm.toymodel(INPUT, PAR)

Su_max = PAR[1]

plt.figure(1)
plt.clf()
plt.subplot(1, 2, 1)
plt.plot(Q[from_index:to_index, pick_catchment], '-b', label='Observed')
plt.plot(QT[from_index:to_index], '-k', label='Simulated')
plt.legend()
plt.ylabel('mm')
plt.xlabel('days')

plt.subplot(1, 2, 2)
plt.plot(Su[from_index:to_index] / Su_max, '.-b', label='Su/Su_max')
plt.ylabel('[-]')
plt.twinx()
plt.plot(Ea[from_index:to_index], '.-r', label='ET')
plt.legend()
plt.ylabel('mm')
plt.xlabel('days')

plt.show()
