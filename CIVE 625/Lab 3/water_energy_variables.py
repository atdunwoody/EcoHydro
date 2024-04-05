def water_energy_variables(option):
    variables = {}
    if option == 1:
        # Forest
        variables['LAI'] = 4.00
        variables['h'] = 20.00
        variables['a'] = 0.12
        variables['g0'] = 15.00
        variables['SMo'] = 80.00
        variables['SMinit'] = 40.00
        variables['Scanmax'] = 4.00
        variables['d'] = 14.644
        variables['z0'] = 1.607
        variables['P'] = 101.20
        variables['zm'] = 22.00
        variables['k'] = 0.4
        variables['esurface'] = 0.95
        variables['s'] = 5.67E-08
        variables['KR'] = 200
        variables['KD1'] = -0.307
        variables['KD2'] = 0.019
        variables['TL'] = 273.00
        variables['T0'] = 293.00
        variables['TH'] = 313.00
        variables['KM1'] = 3.36E-04
        variables['KM2'] = -0.10
        variables['ra'] = 1.23
        variables['cp'] = 1013.00
        variables['gc'] = 1.00
        variables['aT'] = 1.00
    elif option == 2:
        # Grass
        variables['LAI'] = 2.00
        variables['h'] = 0.12
        variables['a'] = 0.23
        variables['g0'] = 30.00
        variables['SMo'] = 40.00
        variables['SMinit'] = 20.00
        variables['Scanmax'] = 2.00
        variables['d'] = 0.077
        variables['z0'] = 0.013
        variables['P'] = 101.2
        variables['zm'] = 2.00
        variables['k'] = 0.4
        variables['esurface'] = 0.95
        variables['s'] = 5.67E-08
        variables['KR'] = 200
        variables['KD1'] = -0.307
        variables['KD2'] = 0.019
        variables['TL'] = 273.00
        variables['T0'] = 293.00
        variables['TH'] = 313.00
        variables['KM1'] = 1.87E-02
        variables['KM2'] = -1.00E-01
        variables['ra'] = 1.23
        variables['cp'] = 1013.00
        variables['gc'] = 1.00
        variables['aT'] = 1.00

    return variables