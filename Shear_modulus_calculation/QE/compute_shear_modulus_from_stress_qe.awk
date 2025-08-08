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
    print "Read:", $1, $5
    stress[$1 + 0] = $5
    Ly[NR] = $10 + 0
}
END {

    printf("Ly = %.6f [Angstrom] \n", Ly[4])

    gamma = (epsilon1 - epsilon3)/Ly[4]
    G1 = -(stress[epsilon1] - stress[epsilon3]) / gamma

    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon1, G1, G1 * conversion_factor)

    gamma = (epsilon2 - epsilon3)/Ly[4]
    G2 = -(stress[epsilon2] - stress[epsilon3]) / gamma

    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon2, G2, G2 * conversion_factor)

    gamma = (epsilon4 - epsilon3)/Ly[4]
    G4 = -(stress[epsilon4] - stress[epsilon3]) / gamma

    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon4, G4, G4 * conversion_factor)

    gamma = (epsilon5 - epsilon3)/Ly[4]
    G5 = -(stress[epsilon5] - stress[epsilon3]) / gamma
    
    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon5, G5, G5 * conversion_factor)
    
    G = (G1 + G2 + G4 + G5)/4
    printf("Average Shear modulus G from energy fit: %.6f Ry/Bohr^3 = %.2f GPa\n", G, G * conversion_factor)
}
