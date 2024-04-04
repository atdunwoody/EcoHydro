def water_energy_variables(option):
    if option == 1:  # Forest
        variables = {
            'LAI': 4.00, 'h': 20.00, 'a': 0.12, 'g0': 15.00, 'SMo': 80.00,
            'SMinit': 40.00, 'Scanmax': 4.00, 'd': 14.644, 'z0': 1.607,
            'P': 101.20, 'zm': 22.00, 'k': 0.4, 'esurface': 0.95, 's': 5.67E-08,
            'KR': 200, 'KD1': -0.307, 'KD2': 0.019, 'TL': 273.00, 'T0': 293.00,
            'TH': 313.00, 'KM1': 3.36E-04, 'KM2': -0.10, 'ra': 1.23,
            'cp': 1013.00, 'gc': 1.00, 'aT': 1.00
        }
    elif option == 2:  # Grass
        variables = {
            'LAI': 2.00, 'h': 0.12, 'a': 0.23, 'g0': 30.00, 'SMo': 40.00,
            'SMinit': 20.00, 'Scanmax': 2.00, 'd': 0.077, 'z0': 0.013,
            'P': 101.2, 'zm': 2.00, 'k': 0.4, 'esurface': 0.95, 's': 5.67E-08,
            'KR': 200, 'KD1': -0.307, 'KD2': 0.019, 'TL': 273.00, 'T0': 293.00,
            'TH': 313.00, 'KM1': 1.87E-02, 'KM2': -1.00E-01, 'ra': 1.23,
            'cp': 1013.00, 'gc': 1.00, 'aT': 1.00
        }

    return variables
