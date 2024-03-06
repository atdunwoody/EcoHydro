import numpy as np

def SF_eq(S0, P, alpha, Tf, SE):
    """
    Fast Flow Reservoir (Sf)
    
    Parameters:
    - S0: Initial storage in the fast flow reservoir
    - P: Precipitation input
    - alpha: Fraction of precipitation that contributes to quick flow
    - Tf: Time parameter for quick flow
    - SE: Saturation excess overland flow contributing to quick flow
    
    Returns:
    - Qf: Quick flow (fast runoff)
    - S_dt: Updated storage in the fast flow reservoir
    - ie: Effective infiltration
    """
    ie = (P * alpha)  # Effective infiltration: part of the precipitation contributing to the fast flow
    S_dt = S0 + (P * alpha) + SE - (S0 / Tf)  # Update storage in the fast flow reservoir with effective infiltration and saturation excess, then subtract outflow
    Qf = S0 / Tf  # Calculate quick flow based on initial storage and time parameter
    return Qf, S_dt, ie

def SU_eq(Su_0, P, Ep, Sumax, alpha, beta):
    """
    Unsaturated Zone Reservoir (Su)
    
    Parameters:
    - Su_0: Initial storage in the unsaturated zone
    - P: Precipitation input
    - Ep: Potential evapotranspiration
    - Sumax: Maximum storage capacity of the unsaturated zone
    - alpha: Fraction of precipitation that does not contribute to immediate runoff
    - beta: Split between recharge and overland flow
    
    Returns:
    - r: Recharge to the saturated zone
    - se: Saturation excess overland flow
    - E: Actual evapotranspiration
    - S_dt: Updated storage in the unsaturated zone
    """
    Su_0 += (P * (1 - alpha))  # Update unsaturated storage with the part of precipitation not contributing to immediate runoff
    
    # Calculate recharge: excess water moves from the unsaturated zone to the saturated zone if storage exceeds capacity
    if Su_0 < Sumax:
        R = 0  # No recharge if storage is below capacity
    else:
        R = Su_0 - Sumax  # Recharge is the excess water above the storage capacity
    
    Su_0 -= R  # Update unsaturated zone storage by subtracting recharge
    
    # Calculate actual evapotranspiration based on current storage and potential evapotranspiration
    E = Ep * (Su_0 / Sumax)
    if E > Ep:
        E = Ep  # Ensure actual evapotranspiration does not exceed potential
    
    Su_0 -= E  # Update unsaturated zone storage by subtracting actual evapotranspiration
    
    S_dt = Su_0  # Final storage in the unsaturated zone
    r = R * beta  # Part of the recharge contributing to the baseflow
    se = R * (1 - beta)  # Part of the recharge contributing to overland flow
    return r, se, E, S_dt

def SS_eq(S0, R, Ts):
    """
    Saturated Zone Reservoir (Ss)
    
    Parameters:
    - S0: Initial storage in the saturated zone
    - R: Recharge from the unsaturated zone
    - Ts: Time parameter for slow flow
    
    Returns:
    - Qs: Slow flow (baseflow)
    - S_dt: Updated storage in the saturated zone
    """
    S_dt = S0 + R - (S0 / Ts)  # Update saturated zone storage with recharge, then subtract outflow
    Qs = S0 / Ts  # Calculate slow flow based on initial storage and time parameter
    return Qs, S_dt

def toymodel(DATA, PAR):
    """
    Toy Model Hydrologic Simulation
    
    Parameters:
    - DATA: Input data matrix with columns [Precipitation, Potential Evapotranspiration]
    - PAR: Model parameters [mir, Su_max, Ts, Tf, beta]
    
    Returns:
    - Ea: Actual evapotranspiration
    - QF: Quick flow (fast runoff)
    - R: Recharge to the saturated zone
    - QS: Slow flow (baseflow)
    - QT: Total flow (QF + QS)
    - Sf: Storage in the fast flow reservoir
    - Su: Storage in the unsaturated zone
    - Ss: Storage in the saturated zone
    - St: Total storage (Su + Ss)
    - AL: Fraction of precipitation that contributes to runoff
    - IE: Effective infiltration
    - SE: Saturation excess overland flow
    """
    M = DATA.shape[0]
    P = DATA[:, 0]  # Precipitation input
    Ep = DATA[:, 1]  # Potential evapotranspiration
    PAR.append(1)  # Add split parameter between R and QFs, set to 1

    # Parameters
    mir = PAR[0]  # Maximum infiltration rate
    Su_max = PAR[1]  # Max storage capacity of unsaturated zone
    Ts = PAR[2]  # Time parameter for slow flow
    Tf = PAR[3]  # Time parameter for quick flow
    beta = PAR[4]  # Split between recharge and overland flow

    # Initialize output arrays
    QF, QS, SE, QT, R, Sf, Su, Ss, Ea, AL, IE = [np.zeros(M) for _ in range(11)]

    # Initial conditions
    S0 = 0.2  # Initial storage for all reservoirs

    # Initial conditions for storage in each reservoir
    Su_dt = S0
    Ss_dt = S0
    Sf_dt = S0 # Assume an initial storage value

    # Main flow routine
    for t in range(M):
        AL[t] = 1 - (1 - np.exp(-P[t] / mir)) / (P[t] / mir) if P[t] > 0 else 0  # Fraction of precipitation contributing to runoff
        
        r, se, E, Su_dt = SU_eq(Su_dt, P[t], Ep[t], Su_max, AL[t], beta)
        R[t], Ea[t], Su[t], SE[t] = r, E, Su_dt, se  # Update recharge, evapotranspiration, unsaturated storage, and saturation excess
        
        Qs, Ss_dt = SS_eq(Ss_dt, r, Ts)
        QS[t], Ss[t] = Qs, Ss_dt  # Update slow flow and saturated zone storage 
    
    for t in range(M):
        Qf, Sf_dt, ie = SF_eq(Sf_dt, P[t], AL[t], Tf, SE[t])
        QF[t], Sf[t], IE[t] = Qf, Sf_dt, ie  # Update quick flow, fast flow reservoir storage, and effective infiltration
    
    QT = QF + QS # Total flow is the sum of quick flow and slow flow
    St = Su + Ss  # Total storage is the sum of unsaturated and saturated zone storages
    
    return Ea, QF, R, QS, QT, Sf, Su, Ss, St, AL, IE, SE

def main():
    r = ''
if __name__ == "__main__":
    main()