# AWK script to compute bulk modulus K from energy vs strain using finite strain method
# Input format: strain energy volume s_xx s_yy s_zz

BEGIN {
    # 1 Ry = 13.605693122994 [eV] * 1.602176634e-19 [J/eV] = 2.179872361e-18 [J]
    # 1 [bohr] = 0.529177210903e-10 [m]
    # 1 [bohr^3] = 1.481847e-31 [m^3]
    # 1 [Ry/Bohr^3] = 2.179872361e-18 / 1.481847e-31 = 1.47105078e13 [Pa] = 1.47105078e4 [GPa]
    conversion_factor = 14710.5  # Ry/Bohr^3 -> GPa
}

{
    print "Read:", $1, $2, $3
    strain[NR] = $1 + 0
    energy[NR] = $2 + 0
    volume[NR] = $3 + 0
}

END {
    # Fit quadratic: E(epsilon) = a * epsilon^2 + b * epsilon + c
    # Use central difference: 
    # (d(dE/d(epsilon))/d(epsilon) is nearly equal to (E(+eps) + E(-eps) - 2E(0)) / (eps^2)
    # K = 1/V0 * (d(dE/d(epsilon))/d(epsilon) = 2a/V0
    
    eps = 0.010
    E_plus = energy[6]
    E_minus = energy[2]
    E_zero = energy[4]
    V0 = volume[4]
    
    second_derivative = (E_plus + E_minus - 2 * E_zero) / (eps^2)
    K = 2*second_derivative / V0
    
    printf("Bulk modulus K from energy fit: %.6f Ry/Bohr^3 = %.2f GPa\n", K, K * conversion_factor)
}
