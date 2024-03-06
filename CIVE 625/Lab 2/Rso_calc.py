import numpy as np

def Rso_calc(Dates, Lat, Lon, Elev):
    """
    Calculate extraterrestrial radiation (Ra) and clear sky solar radiation (Rso).

    Parameters:
    - Dates: numpy array of dates in format [year, month, day].
    - Lat: Latitude in degrees.
    - Lon: Longitude in degrees (not used in calculations directly but included for completeness).
    - Elev: Elevation in meters.

    Returns:
    - Ra: Extraterrestrial radiation in MJ m^-2 day^-1.
    - Rso: Clear sky solar radiation in MJ m^-2 day^-1.
    """

    # Make Julian dates
    years = np.unique(Dates[:, 0])
    jday = np.array([])
    for year in years:
        aux = np.where(Dates[:, 0] == year)[0]
        jday_new = np.arange(1, len(aux) + 1)
        jday = np.concatenate((jday, jday_new))

    # Calculate Ra or So incoming radiation without atmosphere
    dr = 1 + 0.033 * np.cos(2 * np.pi * jday / 365)
    delta = 0.409 * np.sin(2 * np.pi * jday / 365 - 1.39)
    Lat_rad = np.radians(Lat)  # Convert latitude to radians

    Ra = np.zeros((len(jday), len(Lat_rad)))
    Rso = np.zeros((len(jday), len(Lat_rad)))

    for i, lat in enumerate(Lat_rad):
        ws = np.arccos(-np.tan(delta) * np.tan(lat))
        Ra[:, i] = (24 * 60 / np.pi) * 0.08202 * dr * (ws * (np.sin(lat) * np.sin(delta)) + np.sin(ws) * (np.cos(lat) * np.cos(delta)))
        Rso[:, i] = (0.75 + 2 * 10**-5 * Elev[i]) * Ra[:, i]

    return Ra, Rso

def main():
    r =""

if __name__ == "__main__":
    main()
