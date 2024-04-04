def weighted_average_outputs(OUT_A, OUT_B, weight_A, weight_B):
    OUT = {}
    for key in OUT_A.keys():
        if key in ['P', 'Ep']:  # Assuming P and Ep are not to be averaged but taken from OUT_A
            OUT[key] = OUT_A[key]
        else:
            OUT[key] = weight_A * OUT_A[key] + weight_B * OUT_B[key]
    return OUT

