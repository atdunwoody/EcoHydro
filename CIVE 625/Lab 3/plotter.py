import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def Plot_gSM(FC, WP, SMo):
    # Define the function
    def g_SM_calc(FC, WP, SM):
        return np.minimum(np.maximum((SM - WP) / (FC - WP), 0), 1)

    # Define the range of SM values
    SM = np.linspace(0, 2 * SMo, 100)  # Adjust the range as needed
    g_SM = g_SM_calc(FC, WP, SM)
    SM_SMo_ratio = SM / SMo

    # Plotting
    plt.figure(figsize=(10, 5))

    # Plot 1: g_SM vs SM/SMo
    plt.subplot(1, 2, 1)
    plt.plot(SM_SMo_ratio, g_SM, 'LineWidth', 2)
    plt.xlabel('SM / SMo')
    plt.ylabel('g_{SM}')
    plt.title('Plot of g_{SM} versus SM / SMo')
    plt.grid(True)
    plt.axis([0, 2, 0, 1])  # Adjust the axis limits as necessary

    # Plot 2: g_SM vs SM
    plt.subplot(1, 2, 2)
    plt.plot(SM, g_SM, 'LineWidth', 2)
    plt.xlabel('SM')
    plt.ylabel('g_{SM}')
    plt.grid(True)
    plt.axis([0, SMo * 2, 0, 1])  # Adjust the axis limits as necessary

    plt.show()

def plot_FDC(*args):
    plt.figure(4, figsize=(8.27, 1.84))  # Size in inches

    for i, data in enumerate(args):
        QT_sorted = np.sort(data)[::-1]
        permanence = np.arange(1, len(QT_sorted) + 1) / len(QT_sorted)

        if i == 0:
            plt.semilogy(permanence, QT_sorted, '-k', linewidth=1, label='Simulation 1')
        else:
            plt.semilogy(permanence, QT_sorted, '-r', linewidth=1, label='Simulation 2')

    plt.ylabel('Q (mm/d)')
    plt.xlabel('% of Time Exceeded')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()



def plot_FFC(Dates, QT, figureTitle):
    # Ensure QT is a column vector
    QT = np.array(QT).flatten()
    Dates = pd.to_datetime(Dates)
    unique_years = Dates.year.unique()

    annual_max_flows = []

    for year in unique_years:
        streamflow_year = QT[Dates.year == year]
        if len(streamflow_year) > 0:
            annual_max_flows.append(max(streamflow_year))

    sorted_max_flows = np.sort(annual_max_flows)[::-1]
    n_years = len(unique_years)
    exceedance_prob = np.arange(1, n_years + 1) / n_years
    return_periods = 1 / exceedance_prob

    plt.figure()
    plt.loglog(return_periods, sorted_max_flows, 'o-k')
    plt.title(figureTitle)
    plt.xlabel('Return Period (years)')
    plt.ylabel('Annual Maximum Streamflow (mm/d)')
    plt.grid(True)
    plt.show()
