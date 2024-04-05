import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
def plot_gSM(FC, WP, SMo):
    # Define the function
    def g_SM_calc(FC, WP, SM):
        return np.minimum(np.maximum((SM - WP) / (FC - WP), 0), 1)

    # Define the range of SM values
    SM = np.linspace(0, 2 * SMo, 100)  # Adjust the range as needed
    g_SM = g_SM_calc(FC, WP, SM)
    SM_SMo_ratio = SM / SMo

    # Plotting
    plt.figure(figsize=(15, 6))

    # Plot 1: g_SM vs SM/SMo
    plt.subplot(1, 2, 1)
    plt.plot(SM_SMo_ratio, g_SM, linewidth=2)  # Use `linewidth` instead of 'LineWidth'
    plt.xlabel('SM / SMo')
    plt.ylabel('g_{SM}')
    plt.title('Plot of g_{SM} versus SM / SMo')
    plt.grid(True)
    plt.axis([0, 2, 0, 1])  # Adjust the axis limits as necessary

    # Plot 2: g_SM vs SM
    plt.subplot(1, 2, 2)
    plt.plot(SM, g_SM, linewidth=2)  # Use `linewidth` instead of 'LineWidth'
    plt.xlabel('SM')
    plt.ylabel('g_{SM}')
    plt.grid(True)
    plt.axis([0, SMo * 2, 0, 1])  # Adjust the axis limits as necessary

    plt.show()
    
def plot_FFC(Dates, QT, label, ax=None):
    if ax is None:
        fig, ax = plt.subplots()

    # Convert Dates and QT to a pandas Series object
    QT_series = pd.Series(QT, index=Dates)

    # Group by year and get the annual maximum streamflow
    annual_max_flows = QT_series.groupby(QT_series.index.year).max()

    # Sort the annual maximum flows in descending order
    sorted_max_flows = annual_max_flows.sort_values(ascending=False)

    # Calculate exceedance probability and return periods
    n_years = len(sorted_max_flows)
    exceedance_prob = np.arange(1, n_years + 1) / n_years
    return_periods = 1 / exceedance_prob

    # Plotting
    ax.loglog(return_periods, sorted_max_flows, 'o-', label=label)
    ax.set_xlabel('Return Period (years)')
    ax.set_ylabel('Annual Maximum Streamflow (mm/d)')
    ax.grid(True)
    ax.legend()

def plot_FDC(datasets):
    plt.figure(4, figsize=(10, 4))  # Adjust figure size as needed
    plt.clf()
    #unpack label : data from datasets dictionary
    for i, (label, data) in enumerate(datasets.items()):
        QT_sorted = np.sort(data)[::-1]  # Sort data in descending order
        permanence = np.arange(1, len(QT_sorted) + 1) / len(QT_sorted)  # Calculate permanence
        curve_color = '-r' if i == 0 else '-k'  # Set curve color based on index
        plt.semilogy(permanence, QT_sorted, curve_color, linewidth=1, label=label)


    plt.ylabel('Q (mm/d)')
    plt.xlabel('% of Time Exceeded')
    plt.legend()
    plt.grid(True)
    plt.axis('tight')
    plt.show()

