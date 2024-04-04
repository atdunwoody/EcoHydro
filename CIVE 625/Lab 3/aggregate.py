import pandas as pd
import numpy as np

def make_means_new(INPUT, Dates):
    df = pd.DataFrame(INPUT, index=pd.to_datetime(Dates))
    INPUT_MEAN_MONTHLY = df.resample('M').sum() / len(df.index.year.unique())
    INPUT_MEAN_ANNUAL = df.resample('A').sum().mean()
    return INPUT_MEAN_MONTHLY, INPUT_MEAN_ANNUAL

def make_monthly_new(Z, Dates, option):
    df = pd.DataFrame(Z, index=pd.to_datetime(Dates))
    if option == 1:
        Z_MONTHLY = df.resample('M').mean()
    elif option == 2:
        Z_MONTHLY = df.resample('M').sum()
    monthly = Z_MONTHLY.reset_index()
    return monthly, Z_MONTHLY

def make_means(INPUT, Dates):
    df = pd.DataFrame(INPUT, index=pd.to_datetime(Dates))
    INPUT_MEAN_MONTHLY = df.resample('M').mean()
    INPUT_MEAN_ANNUAL = df.resample('A').mean().mean(axis=0)  # Assuming you want the mean across years
    return INPUT_MEAN_MONTHLY, INPUT_MEAN_ANNUAL

def aggregate(OUT, dates):
    date_index = pd.to_datetime(dates)

    # Aggregate Monthly with appropriate option
    monthly_vars = ['Ea', 'QS', 'QF', 'QT', 'IE', 'SE', 'R', 'P', 'Ep', 'Ei', 'Et', 'pot_inf']
    for var in monthly_vars:
        OUT[f"{var}_mon"], _ = make_monthly_new(OUT[var], date_index, 2)

    monthly_means_vars = ['Ss', 'Su', 'St', 'S_canopy']
    for var in monthly_means_vars:
        OUT[f"{var}_mon"], _ = make_monthly_new(OUT[var], date_index, 1)

    # Aggregate mean monthly and mean annual
    mean_vars = ['Ea', 'QS', 'QF', 'QT', 'IE', 'SE', 'R', 'P', 'Ep', 'Ei', 'Et', 'S_canopy', 'pot_inf']
    for var in mean_vars:
        OUT[f"{var}_mean_mon"], OUT[f"{var}_mean_annual"] = make_means_new(OUT[var], date_index)

    # For means that use the original make_means function
    means_vars_direct = ['Ss', 'Su', 'St']
    for var in means_vars_direct:
        OUT[f"{var}_mean_mon"], OUT[f"{var}_mean_annual"] = make_means(OUT[var], date_index)

    return OUT
