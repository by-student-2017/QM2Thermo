
# AWK script to compute shear modulus from shear_results.txt
# Assumes format: strain stress_xy
BEGIN {
    # 1 Ry = 13.605693122994 [eV] * 1.602176634e-19 [J/eV] = 2.179872361e-18 [J]
    # 1 [bohr] = 0.529177210903e-10 [m]
    # 1 [bohr^3] = 1.481847e-31 [m^3]
    # 1 [Ry/Bohr^3] = 2.179872361e-18 / 1.481847e-31 = 1.47105078e13 [Pa] = 1.47105078e4 [GPa]
    #conversion_factor = 14710.5  # Ry/Bohr^3 -> GPa
    conversion_factor = 1.0      # GPa
}
{
    print "Read:", $1, $2, $3, $4, $5, $6, $7, $8, $9
    strain[NR]    = $1 + 0
    energy[NR]    = $2 + 0
    volume[NR]    = $3 + 0
    stress_xx[NR] = $4 + 0
    stress_xy[NR] = $5 + 0
    stress_xz[NR] = $6 + 0
    stress_yy[NR] = $7 + 0
    stress_yz[NR] = $8 + 0
    stress_zz[NR] = $9 + 0
}
END {
    
    # minus
    G1 = (stress_yz[3] - stress_yz[2]) / (2.0*(strain[3] - strain[2]))
    #printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", strain[3], G1, G1 * conversion_factor)
    
    # plus
    G2 = (stress_yz[4] - stress_yz[2]) / (2.0*(strain[4] - strain[2]))
    #printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", strain[4], G2, G2 * conversion_factor)
    
    # minus
    G4 = (stress_yz[5] - stress_yz[2]) / (2.0*(strain[5] - strain[2]))
    #printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", strain[5], G4, G4 * conversion_factor)
    
    # plus
    G5 = (stress_yz[6] - stress_yz[2]) / (2.0*(strain[6] - strain[2]))
    #printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", strain[6], G5, G5 * conversion_factor)
    
    G005 = 0.5*(G1 + G2)
    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", strain[4], G005, G005 * conversion_factor)
    
    G010 = 0.5*(G4 + G5)
    printf("Shear modulus G from energy fit (epislon = %.6f): %.6f Ry/Bohr^3 = %.2f GPa\n", strain[6], G010, G010 * conversion_factor)
    
    G = (G005 + G010)/2
    printf("Average Shear modulus G from energy fit: %.6f Ry/Bohr^3 = %.2f GPa\n", G, G * conversion_factor)
}
