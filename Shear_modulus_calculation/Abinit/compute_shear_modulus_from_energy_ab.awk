# AWK script to compute shear modulus G from energy vs strain using finite strain method
# Input format: strain energy volume stress_xy

BEGIN {
    # 1 Ry = 13.605693122994 [eV] * 1.602176634e-19 [J/eV] = 2.179872361e-18 [J]
    # 1 [bohr] = 0.529177210903e-10 [m]
    # 1 [bohr^3] = 1.481847e-31 [m^3]
    # 1 [Ry/Bohr^3] = 2.179872361e-18 / 1.481847e-31 = 1.47105078e13 [Pa] = 1.47105078e4 [GPa]
    # 1 [eV/Bohr^3] = 1.47105078e4 / 13.605693122994 = 1081.20238102 [GPa]
    conversion_factor = 1081.20238102  # eV/Bohr^3 -> GPa
    
}

{
    strain[NR] = $1 + 0
    energy[NR] = $2 + 0
    volume[NR] = $3 + 0
}

END {
    # Fit quadratic: E(epislon) = a * epsilon^2 + b * epsilon + c
    # Use 5-point fit: -0.01, -0.005, 0.0, 0.005, 0.01
    # Solve for a using central difference approximation:
    # a is nearly equal to (E(+epislon) + E(-epislon) - 2E(0)) / (2*epsilon^2)
    
    #-----------------------------------------------------------------
    eps = 0.010
    E_plus = energy[6]
    E_minus = energy[2]
    E_zero = energy[4]
    V0 = volume[4]
    
    a = (E_plus + E_minus - 2 * E_zero) / (2 * eps^2)
    G = (2 * a) / V0
    #G_GPa = G * conversion_factor
    
    printf("Shear modulus G from energy fit: %.6f eV/Bohr^3 = %.2f GPa\n", G, G * conversion_factor)
    
    #-----------------------------------------------------------------
}
