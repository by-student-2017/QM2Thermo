# AWK script to compute shear modulus from shear_results.txt
# Assumes format: strain stress_xy (=stress_yx)

BEGIN {
    # 1 Ry = 13.605693122994 [eV] * 1.602176634e-19 [J/eV] = 2.179872361e-18 [J]
    # 1 [bohr] = 0.529177210903e-10 [m]
    # 1 [bohr^3] = 1.481847e-31 [m^3]
    # 1 [Ry/Bohr^3] = 2.179872361e-18 / 1.481847e-31 = 1.47105078e13 [Pa] = 1.47105078e4 [GPa]
    #conversion_factor = 14710.5  # Ry/Bohr^3 -> GPa
    conversion_factor = 14710.5  # Ry/Bohr^3 -> GPa

    epsilon0 =  ""
    epsilon1 = -0.010
    epsilon2 = -0.005
    epsilon3 =  0.0
    epsilon4 =  0.005
    epsilon5 =  0.010
}
{
    print "Read:", $1, $8
    stress_yz[$1 + 0] = $8
}
END {

    G1 = -(stress_yz[epsilon1] - stress_yz[epsilon3]) / (epsilon1 - epsilon3)

    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon1, G1, G1 * conversion_factor)

    G2 = -(stress_yz[epsilon2] - stress_yz[epsilon3]) / (epsilon1 - epsilon3)

    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon2, G2, G2 * conversion_factor)

    G4 = -(stress_yz[epsilon4] - stress_yz[epsilon3]) / (epsilon1 - epsilon3)

    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon4, G4, G4 * conversion_factor)

    G5 = -(stress_yz[epsilon5] - stress_yz[epsilon3]) / (epsilon1 - epsilon3)
    
    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon5, G5, G5 * conversion_factor)
    
    G = (G1 + G2 + G4 + G5)/4
    printf("Average Shear modulus G from energy fit: %.6f Ry/Bohr^3 = %.2f GPa\n", G, G * conversion_factor)
}
