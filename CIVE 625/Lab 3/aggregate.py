import numpy as np


def make_means_new(INPUT, Dates):
    years = np.unique(Dates[:, 0])
    years = np.delete(years, 0)  # Remove the first element as in MATLAB

    INPUT_ANNUAL = np.array([np.nansum(INPUT[Dates[:, 0] == year], axis=0) for year in years])
    INPUT_MEAN_ANNUAL = np.nanmean(INPUT_ANNUAL, axis=0)

    INPUT_MEAN_MONTHLY = np.array([np.nansum(INPUT[Dates[:, 1] == month], axis=0) / len(years) for month in range(1, 13)])

    return INPUT_MEAN_MONTHLY, INPUT_MEAN_ANNUAL


def make_means(INPUT, Dates):
    #print(f"Shape of INPUT: {INPUT.shape}")  # Check the shape of INPUT

    years = np.unique(Dates[:, 0])
    years = np.delete(years, 0)
    INPUT_ANNUAL = []

    for year in years:
        aux = np.where(Dates[:, 0] == year)[0]
        #print(f"Indices for year {year}: {aux}")  # Check the indices being used for annual sum
        # Check if INPUT is 1D or 2D
        if INPUT.ndim == 1:
            INPUT_ANNUAL.append(np.nansum(INPUT[aux], axis=0))  # For 1D INPUT
        else:
            INPUT_ANNUAL.append(np.nansum(INPUT[aux, :], axis=0))  # For 2D INPUT

    INPUT_ANNUAL = np.array(INPUT_ANNUAL)
    INPUT_MEAN_ANNUAL = np.nanmean(INPUT_ANNUAL, axis=0)

    INPUT_MEAN_MONTHLY = []
    for i in range(1, 13):  # Adjusted for 0-based indexing
        aux = np.where(Dates[:, 1] == i)[0]
        #print(f"Indices for month {i}: {aux}")  # Check the indices being used for monthly mean
        # Check if INPUT is 1D or 2D
        if INPUT.ndim == 1:
            INPUT_MEAN_MONTHLY.append(np.nanmean(INPUT[aux], axis=0))  # For 1D INPUT
        else:
            INPUT_MEAN_MONTHLY.append(np.nanmean(INPUT[aux, :], axis=0))  # For 2D INPUT

    INPUT_MEAN_MONTHLY = np.array(INPUT_MEAN_MONTHLY)

    return INPUT_MEAN_MONTHLY, INPUT_MEAN_ANNUAL


def make_monthly_new(Z, Dates, option):
    #print(f"Shape of Z: {Z.shape}")  # Confirm the shape of Z

    if len(Z.shape) == 1:
        Z = Z[:, np.newaxis]  # Ensure Z is two-dimensional
    
    years = np.unique(Dates[:, 0])
    #print(f"Unique years: {years}")  # Confirm the unique years

    Z_MONTHLY = []
    for year in years:
        aux1 = np.where(Dates[:, 0] == year)[0]
        #print(f"Number of indices for year {year}: {len(aux1)}")  # Confirm the number of indices per year

        Z_year = Z[aux1, :]
        Dates_year = Dates[aux1]

        Z_monthly = np.empty((12, Z.shape[1]))
        Z_monthly.fill(np.nan)

        for j in range(1, 13):  # Months are from 1 to 12
            aux2 = np.where(Dates_year[:, 1] == j)[0]
            #print(f"Number of indices for month {j} in year {year}: {len(aux2)}")  # Confirm the number of indices per month

            if aux2.size > 0:
                if option == 1:
                    Z_monthly[j - 1] = np.mean(Z_year[aux2, :], axis=0)
                elif option == 2:
                    Z_monthly[j - 1] = np.sum(Z_year[aux2, :], axis=0)

        Z_MONTHLY.append(Z_monthly)

    Z_MONTHLY = np.vstack(Z_MONTHLY)
    return Z_MONTHLY


def aggregate(OUT, dates):
    # Aggregate Monthly
    OUT['E_mon'] = make_monthly_new(OUT['Ea'], dates, 2)
    OUT['Ss_mon'] = make_monthly_new(OUT['Ss'], dates, 1)
    OUT['Su_mon'] = make_monthly_new(OUT['Su'], dates, 1)
    OUT['St_mon'] = make_monthly_new(OUT['St'], dates, 1)
    OUT['QS_mon'] = make_monthly_new(OUT['QS'], dates, 2)
    OUT['QF_mon'] = make_monthly_new(OUT['QF'], dates, 2)
    OUT['QT_mon'] = make_monthly_new(OUT['QT'], dates, 2)
    OUT['IE_mon'] = make_monthly_new(OUT['IE'], dates, 2)
    OUT['SE_mon'] = make_monthly_new(OUT['SE'], dates, 2)
    OUT['R_mon'] = make_monthly_new(OUT['R'], dates, 2)
    OUT['P_mon'] = make_monthly_new(OUT['P'], dates, 2)
    OUT['PET_mon'] = make_monthly_new(OUT['Ep'], dates, 2)
    OUT['Ei_mon'] = make_monthly_new(OUT['Ei'], dates, 2)
    OUT['Et_mon'] = make_monthly_new(OUT['Et'], dates, 2)
    OUT['S_canopy_mon'] = make_monthly_new(OUT['S_canopy'], dates, 1)
    OUT['pot_inf_mon'] = make_monthly_new(OUT['pot_inf'], dates, 2)

    # Aggregate mean_monthly and mean annual
    OUT['E_mean_mon'], OUT['E_mean_annual'] = make_means_new(OUT['Ea'], dates)
    OUT['Ss_mean_mon'], OUT['Ss_mean_annual'] = make_means(OUT['Ss'], dates)
    OUT['Su_mean_mon'], OUT['Su_mean_annual'] = make_means(OUT['Su'], dates)
    OUT['St_mean_mon'], OUT['St_mean_annual'] = make_means(OUT['St'], dates)
    OUT['QS_mean_mon'], OUT['QS_mean_annual'] = make_means_new(OUT['QS'], dates)
    OUT['QF_mean_mon'], OUT['QF_mean_annual'] = make_means_new(OUT['QF'], dates)
    OUT['QT_mean_mon'], OUT['QT_mean_annual'] = make_means_new(OUT['QT'], dates)
    OUT['IE_mean_mon'], OUT['IE_mean_annual'] = make_means_new(OUT['IE'], dates)
    OUT['SE_mean_mon'], OUT['SE_mean_annual'] = make_means_new(OUT['SE'], dates)
    OUT['R_mean_mon'], OUT['R_mean_annual'] = make_means_new(OUT['R'], dates)
    OUT['P_mean_mon'], OUT['P_mean_annual'] = make_means_new(OUT['P'], dates)
    OUT['PET_mean_mon'], OUT['PET_mean_annual'] = make_means_new(OUT['Ep'], dates)
    OUT['Ei_mean_mon'], OUT['Ei_mean_annual'] = make_means_new(OUT['Ei'], dates)
    OUT['Et_mean_mon'], OUT['Et_mean_annual'] = make_means_new(OUT['Et'], dates)
    OUT['S_canopy_mean_mon'], OUT['S_canopy_mean_annual'] = make_means_new(OUT['S_canopy'], dates)
    OUT['pot_inf_mean_mon'], OUT['pot_inf_mean_annual'] = make_means_new(OUT['pot_inf'], dates)

    return OUT
