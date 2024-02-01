import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt

class AirParcel:
    def __init__(self, temperature, pressure =101.3, RH = None, psy = None):
        self.T = temperature  # Temperature in degrees Celsius
        self.P = pressure  # Pressure in kPa
        self.RH = RH  # Relative humidity as a fraction
        self.c_p = 1.013  # kJ/kg/K, specific heat of dry air at constant pressure
        self.R_v = None  #0.4615 kJ/kg/K, specific gas constant for water vapor
        self.R_d = None  #0.287 kJ/kg/K, specific gas constant for dry air
        self.rho_v = None # Vapor density in g/m^3
        self.rho_d = None # Dry air density in g/m^3
        self.e_sat = self.calc_esat() # Saturation vapor pressure in kPa
        self.e = self.calc_e()  # Vapor pressure in kPa
        self.LH = self.calc_LH() # Latent heat of vaporization in MJ/kg
        self.psy = psy if psy is not None else self.calc_psy() # Psychrometric constant in kPa/°C
        self.T_dew = self.calc_dewpoint() # Dewpoint temperature in degrees Celsius
        self.T_wet = self.calc_wetbulb() # Wetbulb temperature in degrees Celsius
        self.specific_humidity = self.calc_specific_humidity() # Specific humidity in g/g
        self.mixing_ratio = None # Mixing ratio in g/g
       
    def __repr__(self):
        return f"AirParcel({self.T}, {self.P})"

    def calc_rho_v(self):
        """
        Calculate vapor density from pressure, specific gas constant for water vapor, and temperature
        P: pressure in kPa
        R_v: specific gas constant for water vapor in kJ/kg/K
        T: temperature in degrees Celsius
        return: vapor density in g/m^3
        """
        rho_v = self.P / (self.R_v * (self.T + 273.15))
        return rho_v
    
    def calc_rho_d(self):
        """
        Calculate dry air density from pressure, specific gas constant for dry air, and temperature
        P: pressure in kPa
        R_d: specific gas constant for dry air in kJ/kg/K
        T: temperature in degrees Celsius
        return: dry air density in kg/m^3
        """
        rho_d = self.P / (self.R_d * (self.T + 273.15))
        return rho_d

    def calc_e(self):
        """
        Calculate vapor pressure from temperature
        Uses instance attribute for temperature (T)
        """
        if self.RH is not None:
            e = self.RH * self.e_sat
        elif self.rho_v is not None:
            e = self.rho_v*self.R_v*self.T
        else:
            print("Need to provide RH or vapor density to calculate vapor pressure")
            return None
        return e    

    def calc_esat(self):
        """
        Calculate saturation vapor pressure from temperature
        T: temperature in degrees Celsius
        return: saturation vapor pressure in kPa
        """
        esat = 0.6108 * np.exp(17.27 * self.T / (self.T + 237.3))
        return esat

    def calc_LH(self):
        """
        LH = Lambda #print("\u03BB")
        Calculate latent heat of vaporization from temperature
        T: temperature in degrees Celsius
        return: latent heat of vaporization (lambda) in MJ/kg
        """
        LH = 2.501 - 0.002361 * self.T
        return LH

    def calc_RH(self):
        """
        Calculate relative humidity from vapor pressure and saturation vapor pressure
        e: vapor pressure in kPa
        e_sat: saturation vapor pressure in kPa
        return: relative humidity in percent
        """
        RH = (self.e / self.e_sat)
        return RH
    
    def calc_psy(self):
        """
        γ = psychrometric constant
        Calculate psychrometric constant from temperature
        T: temperature in degrees Celsius
        """
        psychrometric_constant = self.c_p * self.P / (0.622 * self.calc_LH(self.T))
        return psychrometric_constant

    def calc_dewpoint(self):
        """
        Calculate dewpoint temperature (T_dew) from vapor pressure
        e: vapor pressure in kPa
        return: dewpoint temperature in degrees Celsius
        """
        T_dew = (np.log(self.e) + 0.49299) / (0.0707 - 0.00421 * np.log(self.e))
        return T_dew

    def calc_wetbulb(self):
        """
        Calculate wetbulb temperature (T_wet) from vapor pressure, psychrometric constant, dry bulb temperature, and saturation vapor pressure
        e: vapor pressure in kPa
        T_dry: dry bulb temperature in degrees Celsius
        e_sat: saturation vapor pressure in kPa
        return: wetbulb temperature in degrees Celsius
        """
        T_dry = self.T # Dry bulb temperature is equal to air temperature
        T_wet = (self.e + self.psy * T_dry) / (self.e_sat + self.psy)
        return T_wet

    def calc_specific_humidity(self):
        """
        specific humidity = q
        Calculate specific humidity from vapor density and dry air density in moist air sample
        rho_v: vapor density in g/m^3 
        rho_d: dry air density in kg/m^3
        return: specific humidity in g/kg
        """
        if self.rho_v is not None and self.rho_d is not None:    
            q = self.rho_v / (self.rho_v + self.rho_d)
        else:
            q = 0.622 * self.e / (self.P - 0.38 * self.e)
        return q

    def calc_mixing_ratio(self):
        """
        Mixing ratio = r
        Calculate mixing ratio from vapor density and dry air density in moist air sample
        rho_v: vapor density in g/m^3 
        rho_d: dry air density in kg/m^3
        return: mixing ratio in g/kg
        """
        r = self.rho_v / self.rho_d
        return r

class SolarRadiation:
    def __init__(self, latitude, day_of_year, time_of_day = None, solar_constant = 1367):
        self.latitude = latitude # Latitude in degrees
        self.day_of_year = day_of_year  # Day of year
        self.time_of_day = time_of_day   # Time of day in hours
        self.solar_constant = solar_constant # Solar constant in W/m^2
        self.d_r = self.calc_eccentricity_factor() # Eccentricity factor
        self.sigma = self.calc_solar_declination() # Solar declination in degrees
        self.omega = None # Hour angle in degrees
        self.omega_s = self.calc_sunset_hour_angle() # Sunset hour angle in degrees

    def __repr__(self):
        return f"Radiation({self.latitude}, {self.day_of_year})"
    
    def calc_eccentricity_factor(self):
        """
        Calculate eccentricity factor from day of year
        day_of_year: day of year
        return: eccentricity factor
        """
        d_r = 1 + 0.033 * np.cos(2 * np.pi * self.day_of_year / 365)
        return d_r

    def calc_solar_declination(self):
        """
        Calculate solar declination from day of year
        day_of_year: day of year
        return: solar declination in degrees
        """
        sigma= 0.4093 * np.sin(2 * np.pi * self.day_of_year / 365 - 1.405)
        return sigma
    
    def calc_hour_angle(self):
        """
        Calculate hour angle from time of day
        time_of_day: time of day in hours
        return: hour angle in degrees
        """
        omega = np.pi* (12 - self.time_of_day) / 12
        return omega

    def calc_sunset_hour_angle(self):
        """
        Calculate sunset hour angle from latitude and solar declination
        latitude: latitude in degrees
        sigma: solar declination in degrees
        return: sunset hour angle in degrees
        """
        sigma = self.calc_solar_declination()
        omega_s = np.arccos(-np.tan(np.radians(self.latitude)) * np.tan(sigma))
        return omega_s

    def radiation_flux(self):
        """
        Calculate radiation flux from latitude, day of year, and time of day
        latitude: latitude in degrees
        day_of_year: day of year
        time_of_day: time of day in hours
        return: radiation flux in W / m^2
        """
        Sd_o = (37.7 * self.d_r * 
        (self.omega_s * np.sin(np.radians(self.latitude)) * np.sin(self.sigma) + 
            np.cos(np.radians(self.latitude)) * np.cos(self.sigma) * np.sin(self.omega_s)))
        sec_per_day = 24*60*60
        radiation_flux = Sd_o / sec_per_day * 10**6 # Convert from MJ/m^2/day to W/m^2
        return radiation_flux

class SurfaceRadiation:
    def __init__(self): 
        self.RH = None # Relative humidity as a fraction
        self.T_a = None # Air temperature in degrees Celsius
        self.T_s = None # Surface temperature in degrees Celsius
        self.emmisivity = None # Emmisivity of the surface
        self.S_in = None # Incident short wave solar radiation in cal / cm^2 / day
        self.LW = None # Long wave radiation in W/m^2
        self.R_n = None # Net radiation in W/m^2
        self.albedo = None # Albedo of the surface
        self.G = None # Heat flux into the ground in W/m^2
        self.H = None # Turbulent sensible heat flux in W/m^2
        self.E = None # Evaporation rate in m/s
        self.sigma = 5.67 * (10**-8) #stefan boltzmann constant

    def calc_LW(self):
        """
        Calculate long wave radiation from surface temperature and emmisivity
        T_s: surface temperature in degrees Celsius
        emmisivity: emmisivity of the surface
        return: long wave radiation in W/m^2
        """
        LW = self.emmisivity * self.sigma * (self.T_s + 273.15)**4
        return LW
    
    def calc_R_n(self):
        """
        Calculate net radiation from incident short wave solar radiation and relative humidity
        S_in: incident short wave solar radiation in W/m^2
        return: net radiation in W/m^2
        """
        R_n = self.S_in * (1 - self.albedo) - self.LW
        return R_n
    
    def calc_G(self):
        """
        Calculate heat flux into the ground from net radiation
        R_n: net radiation in W/m^2
        return: heat flux into the ground in W/m^2
        """
        G = self.R_n - self.H - self.E
        return G
    
    def calc_E(self):
        """
        Calculate evaporation rate from net radiation
        R_n: net radiation in W/m^2
        return: evaporation rate in m/s
        """
        E = self.R_n - self.H - self.G
        return E
    
    def calc_H(self):
        """
        Calculate turbulent sensible heat flux from net radiation
        R_n: net radiation in W/m^2
        return: turbulent sensible heat flux in W/m^2
        """
        H = self.R_n - self.E - self.G
        return H
#============================Question 1=======================================

parcel = AirParcel(32, 95, RH = 0.5, psy = 0.066)  # 32 degrees Celsius, 95 kPa, 50% RH, 0.066 kPa/°C

T_dew = parcel.calc_dewpoint()
T_wet = parcel.calc_wetbulb()
specific_humidity = parcel.calc_specific_humidity()

#Format to 2 decimal places
print("Question 1")
print(f"a. Dewpoint Temp {T_dew:.2f}°C")
print(f"b. Wetbulb Temp: {T_wet:.2f}°C")
print(f"c. Specific Humidity: {specific_humidity:.4f} g/kg\n")

#============================Question 3=======================================

fc = 40.585258 # Latitude of Fort Collins, CO
ep = 31.761877 # Latitude of El Paso, TX
day = 38 #February 7th is the 38th day of the year

fc_rad = SolarRadiation(fc, day)  # 40 degrees latitude, on February 7th
ep_rad = SolarRadiation(ep, day)  # 31 degrees latitude, on February 7th
fc_rad.latitude = 100

fc_yearly_rad =[]
ep_yearly_rad = []
#Create a list of daily radiation flux for Fort Collins and El Paso
for i in range(365):
    fc_rad = SolarRadiation(fc, i)
    ep_rad = SolarRadiation(ep, i)
    fc_yearly_rad.append(fc_rad.radiation_flux())
    ep_yearly_rad.append(ep_rad.radiation_flux())

#Plot the radiation flux for Fort Collins and El Paso
sns.set()  # This sets the seaborn style
plt.figure(figsize=(10, 6))  # You can adjust the figure size as needed

# Creating lineplots
sns.lineplot(x=range(365), y=fc_yearly_rad, label="Fort Collins")
sns.lineplot(x=range(365), y=ep_yearly_rad, label="El Paso")

# Adding title and labels
plt.title("Daily Extraterrestrial Solar Radiation ")
plt.xlabel("Day of Year")
plt.ylabel(r"Incident Radiation (W/m$^2$)")

# Setting the axis limits
plt.xlim(0, 365)
plt.ylim(0, 1000)

# Showing the legend and plot
plt.legend()
plt.show()

#============================Question 4=======================================

surf = SurfaceRadiation()

day_to_sec = 24*60*60 # Conversion factor from days to seconds
cal_to_W = 4.184 # Conversion factor from cal to W
cm2_to_m2 = 10**-4 # Conversion factor from cm^2 to m^2
S_in = 468 # Incident short wave solar radiation in cal / cm^2 / day
surf.S_in = 468 * cal_to_W * cm2_to_m2 / day_to_sec # Conversion factor from cal/cm^2/day to W/m^2
surf.RH = 0.66 # Relative humidity as a fraction
surf.T_a = 17.94 # Air temperature in degrees Celsius
surf.T_s = surf.T_a # Surface temperature is equal to air temperature
surf.emmisivity = 0.9 # From Table 5.2 used lower end of range for short grass
surf.albedo = 0.23 # From Table 5.2 Albedo for short grass chosen from middle of range
surf.LW = surf.calc_LW()
surf.R_n = surf.calc_R_n()
print("Question 4")
print(f"Net Radiation: {surf.R_n:.0f} W/m^2\n")

#============================Question 5=======================================
"""
Consider the following observations: net radiation Rn = 200 Wm-2; heat flux into the ground, G = 40 Wm-2; and evaporation rate, E = 5 10-8 ms-1.
a)	Calculate the turbulent sensible heat flux, H, in Wm-2
b)	Was the atmosphere stable or unstable? Why?
c)	Was the soil warming up or cooling down? Why?

"""

E = 5 * 10**-8 # Evaporation rate in m/s
E = E * 1000 / (24 * 60 * 60)  # Conversion factor from m/s to mm/day
E = E * 28.6 # 1 mm/day is equivalent to 28.6 W m−2
surf.E = E
surf.R_n = 200 # Net radiation in W/m^2
surf.G = 40 # Heat flux into the ground in W/m^2
surf.H = surf.calc_H()
#Answer to a
print("Question 5")
print(f"Turbulent Sensible Heat Flux: {surf.H:.2f} W/m^2")

#Answer to b
#The atmosphere is unstable because the net radiation is positive, which means the surface is warmer than the air above it
# This will cause the air to rise and the atmosphere to be unstable.

#Answer to c
#The soil is warming up because the net radiation is positive and the heat flux into the ground is positive.