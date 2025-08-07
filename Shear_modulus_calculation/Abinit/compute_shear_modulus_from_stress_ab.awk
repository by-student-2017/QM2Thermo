
# AWK script to compute shear modulus from shear_results.txt
# Assumes format: strain stress_xy
BEGIN {
    # 1 Ry = 13.605693122994 [eV] * 1.602176634e-19 [J/eV] = 2.179872361e-18 [J]
    # 1 [bohr] = 0.529177210903e-10 [m]
    # 1 [bohr^3] = 1.481847e-31 [m^3]
    # 1 [Ry/Bohr^3] = 2.179872361e-18 / 1.481847e-31 = 1.47105078e13 [Pa] = 1.47105078e4 [GPa]
    #conversion_factor = 14710.5  # Ry/Bohr^3 -> GPa
    conversion_factor = 1.0      # GPa
    
    epsilon0 =  ""
    epsilon1 = -0.01
    epsilon2 = -0.005
    epsilon3 =  0.0
    epsilon4 =  0.005
    epsilon5 =  0.01
}
{
    print "Read:", $1, $5
    strain[$1 + 0] = $5
}
END {
    # Central difference approximation: G = (sigma(+e) - sigma(-e)) / (2 * e)
    if ((strain[epsilon5] != "") && (strain[epsilon1] != "")) {
        G1 = (strain[epsilon5] - strain[epsilon1]) / (2 * epsilon5)
        printf("Shear modulus (G) from +/- %.3f strain: %.6f GPa\n", epsilon5, G1*conversion_factor)
    }
    if ((strain[epsilon4] != "") && (strain[epsilon2] != "")) {
        G2 = (strain[epsilon4] - strain[epsilon2]) / (2 * epsilon4)
        printf("Shear modulus (G) from +/- %.3f strain: %.6f GPa\n", epsilon4, G2*conversion_factor)
    }
    printf("Average Shear modulus (G): %.6f GPa\n", (G1+G2)/2*conversion_factor)
}
