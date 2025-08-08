# AWK script to compute bulk modulus from bulk_results.txt
# Assumes format: strain energy volume s_xx s_yy s_zz

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
    print "Read:", $1, $4, $7, $9
    pressure[$1 + 0] = - ($4 + $7 + $9) / 3.0  # -trace(s)/3
    volume[NR] = $3 + 0
}
END {

    V0 = volume[4]

    kappa = (volume[2] - V0)/V0
    K1 = (pressure[epsilon1] - pressure[epsilon3]) / kappa

    printf("Bulk modulus K from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon1, K1, K1 * conversion_factor)

    kappa = (volume[3] - V0)/V0
    K2 = (pressure[epsilon2] - pressure[epsilon3]) / kappa

    printf("Bulk modulus K from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon2, K2, K2 * conversion_factor)

    kappa = (volume[5] - V0)/V0
    K4 = (pressure[epsilon4] - pressure[epsilon3]) / kappa

    printf("Bulk modulus K from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon4, K4, K4 * conversion_factor)

    kappa = (volume[6] - V0)/V0
    K5 = (pressure[epsilon5] - pressure[epsilon3]) / kappa
    
    printf("Bulk modulus K from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", epsilon5, K5, K5 * conversion_factor)
    
    K = (K1 + K2 + K4 + K5)/4
    printf("Average Bulk modulus K from energy fit: %.6f Ry/Bohr^3 = %.2f GPa\n", K, K * conversion_factor)
}
