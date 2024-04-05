def weighted_average_outputs(OUT_A, OUT_B, weight_A, weight_B):
    OUT = {}
    fields = ['Ea', 'QF', 'R', 'QS', 'QT', 'Sf', 'Su', 'Ss', 'St', 'AL', 'IE', 'SE', 'Ei', 'Et', 'S_canopy', 'pot_inf']
    for field in fields:
        OUT[field] = weight_A * OUT_A[field] + weight_B * OUT_B[field]

    # Assuming these are constants across A and B
    OUT['P'] = OUT_A['P']
    OUT['Ep'] = OUT_A['Ep']

    return OUT
