import sympy as sp

# Define the symbols
e_sat, T = sp.symbols('e_sat T')

# Define the equation
equation = sp.Eq(0.6108 * sp.exp((17.27 * T) / (T + 237.3)), e_sat)

# Solve the equation for T
solution = sp.solve(equation, T, dict=True)
print(solution)