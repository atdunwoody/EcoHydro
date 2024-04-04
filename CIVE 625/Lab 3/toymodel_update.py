import numpy as np


def SU_eq_mod(Su_0, Precip, Sumax, alpha, beta, vars, Sin, Lnet, Temp, u2, e_actual, S_canopy_old):
    def e_sat_calc(Temp):
        return 0.6108 * np.exp((17.27 * Temp) / (237.3 + Temp))

    def ra_calc(u, zm, d, z0, zm_p):
        return 1 / (0.4**2 * u) * np.log((zm - d) / z0) * np.log((zm_p - d) / (z0 / 10))

    def g_r_calc(Sin, KR):
        return max(Sin * (1000 + KR) / (1000 * (Sin + KR)), 0)

    def g_D_calc(KD1, KD2, VPD):
        return max(1 + KD1 * VPD + KD2 * VPD**2, 0)

    def g_T_calc(Temp, TL, TH, T0, alpha_T):
        return max((Temp - TL) * (TH - Temp)**alpha_T / ((T0 - TL) * (TH - T0)**alpha_T), 0)

    def g_SM_calc(FC, WP, SM):
        return min(max((SM - WP) / (FC - WP), 0), 1)

    zm = vars['zm']
    z0 = vars['z0']
    d = vars['d']
    g_0 = vars['g0']
    KR = vars['KR']
    KD1 = vars['KD1']
    KD2 = vars['KD2']
    TH = vars['TH']
    T0 = vars['T0']
    TL = vars['TL']
    alpha_gt = (TH - T0) / (T0 - TL)
    FC = vars['FC']
    WP = vars['WP']
    albedo = vars['a']
    cp = vars['cp']
    rho_air = vars['ra']
    S_canopy_max = vars['Scanmax']
    Patm = vars['Patm']

    S_canopy_temp = S_canopy_old + Precip
    rain_pass = S_canopy_temp - S_canopy_max if S_canopy_temp > S_canopy_max else 0
    S_canopy_temp = min(S_canopy_temp, S_canopy_max)

    Su_0 += rain_pass * (1 - alpha)
    R = Su_0 - Sumax if Su_0 >= Sumax else 0
    Su_0 -= R

    r_a = ra_calc(u2, zm, d, z0, zm)
    e_sat = e_sat_calc(Temp)
    VPD = e_sat - e_actual
    Temp_K = Temp + 273.15

    g_r = g_r_calc(Sin, KR)
    g_D = g_D_calc(KD1, KD2, VPD)
    g_T = g_T_calc(Temp_K, TL, TH, T0, alpha_gt)
    g_SM = g_SM_calc(FC, WP, Su_0)

    g_s = g_0 * g_r * g_D * g_T * g_SM
    r_s = 1000 / g_s if g_s != 0 else 1e6

    Rn = Sin * (1 - albedo) + Lnet
    Lambda = 1000 * (2.501 - 0.002361 * Temp)
    Delta = 4098 * e_sat / (237.3 + Temp)**2
    Gamma = cp * Patm / (0.622 * Lambda * 1000)
    Wm2_mm = 1 / (Lambda * 1000) * 86400

    A = Rn
    E_pot = (Delta * A + rho_air * cp * VPD / r_a) / (Delta + Gamma * (1 + r_s / r_a))
    E_pot_mm = E_pot * Wm2_mm
    excess_ET = E_pot_mm - S_canopy_temp if E_pot_mm > S_canopy_temp else 0
    E_int_mm = min(E_pot_mm, S_canopy_temp)
    E_tran_mm = (E_pot - E_int_mm) * Wm2_mm

    S_canopy = S_canopy_temp - E_int_mm
    Su_0 -= E_tran_mm
    E_total_mm = E_int_mm + E_tran_mm

    S_dt = Su_0
    r = R * beta
    se = R * (1 - beta)

    return r, se, S_dt, S_canopy, E_int_mm, E_tran_mm, E_total_mm, rain_pass, r_a, r_s

def SF_eq(S0, P, alpha, Tf, SE):
    ie = P * alpha
    S_dt = S0 + ie + SE - (S0 / Tf)
    Qf = S0 / Tf
    return Qf, S_dt, ie

def SS_eq(S0, R, Ts):
    S_dt = S0 + R - (S0 / Ts)
    Qs = S0 / Ts
    return Qs, S_dt

def toymodel_update(DATA, PAR, vars):
    M = len(DATA)
    P = DATA[:, 0]
    Ep = DATA[:, 1]
    e = DATA[:, 2]
    u2 = DATA[:, 3]
    Sin = DATA[:, 4]
    Temp = DATA[:, 5]
    Lnet = DATA[:, 6]

    mir = PAR[0]
    Su_max = PAR[1]
    Ts = PAR[2]
    Tf = PAR[3]
    beta = 1

    QF = np.zeros(M)
    QS = np.zeros(M)
    SE = np.zeros(M)
    QT = np.zeros(M)
    R = np.zeros(M)
    Sf = np.zeros(M)
    Su = np.zeros(M)
    Ss = np.zeros(M)
    Ea = np.zeros(M)
    AL = np.zeros(M)
    IE = np.zeros(M)
    Ei = np.zeros(M)
    Et = np.zeros(M)
    S_canopy = np.zeros(M)
    pot_inf = np.zeros(M)

    S0 = 200
    Su_dt = S0
    Ss_dt = S0
    Sf_dt = S0
    S_canopy_old = 0

    ET_vars = {'ra': np.zeros(M), 'rs': np.zeros(M)}

    for t in range(M):
        if P[t] > 0:
            AL[t] = 1 - (1 - np.exp(-P[t] / mir)) / (P[t] / mir)
        else:
            AL[t] = 0

        r, se, Su_dt, S_can, E_int_mm, E_tran_mm, E_total_mm, rain_pass, ra, rs = SU_eq_mod(Su_dt, P[t], Su_max, AL[t], beta, vars, Sin[t], Lnet[t], Temp[t], u2[t], e[t], S_canopy_old)

        ET_vars['ra'][t] = ra
        ET_vars['rs'][t] = rs

        R[t] = r
        Ea[t] = E_total_mm
        Su[t] = Su_dt
        SE[t] = se
        Ei[t] = E_int_mm
        Et[t] = E_tran_mm
        S_canopy[t] = S_can
        S_canopy_old = S_can
        pot_inf[t] = rain_pass

        Qs, Ss_dt = SS_eq(Ss_dt, r, Ts)

        QS[t] = Qs
        Ss[t] = Ss_dt

    for t in range(M):
        Qf, Sf_dt, ie = SF_eq(Sf_dt, P[t], AL[t], Tf, SE[t])

        QF[t] = Qf
        Sf[t] = Sf_dt
        IE[t] = ie

    QT = QF + QS
    St = Su + Ss

    OUT = {
        'Ea': Ea, 'QF': QF, 'R': R, 'QS': QS, 'QT': QT,
        'Sf': Sf, 'Su': Su, 'Ss': Ss, 'St': St, 'AL': AL, 'IE': IE,
        'SE': SE, 'Ei': Ei, 'Et': Et, 'S_canopy': S_canopy,
        'pot_inf': pot_inf, 'ET_vars': ET_vars, 'P': P, 'Ep': Ep
    }

    return OUT