# AWK script to compute shear modulus from elastic_results.txt
# Assumes format: strain stress_xx, ...

BEGIN {
    # 1 Ry = 13.605693122994 [eV] * 1.602176634e-19 [J/eV] = 2.179872361e-18 [J]
    # 1 [bohr] = 0.529177210903e-10 [m]
    # 1 [bohr^3] = 1.481847e-31 [m^3]
    # 1 [Ry/Bohr^3] = 2.179872361e-18 / 1.481847e-31 = 1.47105078e13 [Pa] = 1.47105078e4 [GPa]
    #conversion_factor = 14710.5  # Ry/Bohr^3 -> GPa
    conversion_factor = 14710.5  # Ry/Bohr^3 -> GPa
}
{
    print "Read:", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
    strain[NR]    = $1 + 0
    energy[NR]    = $2 + 0
    volume[NR]    = $3 + 0
    stress_xx[NR] = $4 + 0
    stress_xy[NR] = $5 + 0
    stress_xz[NR] = $6 + 0
    stress_yy[NR] = $7 + 0
    stress_yz[NR] = $8 + 0
    stress_zz[NR] = $9 + 0
    Lx[NR]        = $10 + 0
    Ly[NR]        = $11 + 0
    Lz[NR]        = $12 + 0
}
END {
    # elastic_results.txt
    # Item                d1(s_xx), ...
    # Standard: +0.0000
    # dir 1 minus
    # dir 1 plus
    # dir 2 minus
    # dir 2 plus
    # dir 3 minus
    # ...
    
    
    for (i=1; i<=6; i++) {
      
      #printf("\n strain: %15.8f \n", strain[3+2*(i-1)])
      # minues
      C1neg = -(stress_xx[3+2*(i-1)] - stress_xx[2]) / (strain[3+2*(i-1)] / Lx[2]) * conversion_factor # d1
      C2neg = -(stress_yy[3+2*(i-1)] - stress_yy[2]) / (strain[3+2*(i-1)] / Ly[2]) * conversion_factor # d2
      C3neg = -(stress_zz[3+2*(i-1)] - stress_zz[2]) / (strain[3+2*(i-1)] / Lz[2]) * conversion_factor # d3
      C4neg = -(stress_yz[3+2*(i-1)] - stress_yz[2]) / (strain[3+2*(i-1)] / Lz[2]) * conversion_factor # d4
      C5neg = -(stress_xz[3+2*(i-1)] - stress_xz[2]) / (strain[3+2*(i-1)] / Lz[2]) * conversion_factor # d5
      C6neg = -(stress_xy[3+2*(i-1)] - stress_xy[2]) / (strain[3+2*(i-1)] / Ly[2]) * conversion_factor # d6
      
      #printf("\n strain: %15.8f \n", strain[4+2*(i-1)])
      # plus
      C1pos = -(stress_xx[4+2*(i-1)] - stress_xx[2]) / (strain[4+2*(i-1)] / Lx[2]) * conversion_factor # d1
      C2pos = -(stress_yy[4+2*(i-1)] - stress_yy[2]) / (strain[4+2*(i-1)] / Ly[2]) * conversion_factor # d2
      C3pos = -(stress_zz[4+2*(i-1)] - stress_zz[2]) / (strain[4+2*(i-1)] / Lz[2]) * conversion_factor # d3
      C4pos = -(stress_yz[4+2*(i-1)] - stress_yz[2]) / (strain[4+2*(i-1)] / Lz[2]) * conversion_factor # d4
      C5pos = -(stress_xz[4+2*(i-1)] - stress_xz[2]) / (strain[4+2*(i-1)] / Lz[2]) * conversion_factor # d5
      C6pos = -(stress_xy[4+2*(i-1)] - stress_xy[2]) / (strain[4+2*(i-1)] / Ly[2]) * conversion_factor # d6
      
      # dir = 1, C1${dir}...C6${dir}
      C1[i] = 0.5 * (C1neg + C1pos)
      C2[i] = 0.5 * (C2neg + C2pos)
      C3[i] = 0.5 * (C3neg + C3pos)
      C4[i] = 0.5 * (C4neg + C4pos)
      C5[i] = 0.5 * (C5neg + C5pos)
      C6[i] = 0.5 * (C6neg + C6pos)
    }
    
    # Output final values
    
    C11all = C1[1]
    C22all = C2[2]
    C33all = C3[3]
    
    C12all = 0.5*(C1[2]+C2[1])
    C13all = 0.5*(C1[3]+C3[1])
    C23all = 0.5*(C2[3]+C3[2])
    
    C44all = C4[4]
    C55all = C5[5]
    C66all = C6[6]
    
    C14all = 0.5*(C1[4]+C4[1])
    C15all = 0.5*(C1[5]+C5[1])
    C16all = 0.5*(C1[6]+C6[1])
    
    C24all = 0.5*(C2[4]+C4[2])
    C25all = 0.5*(C2[5]+C5[2])
    C26all = 0.5*(C2[6]+C6[2])
    
    C34all = 0.5*(C3[4]+C4[3])
    C35all = 0.5*(C3[5]+C5[3])
    C36all = 0.5*(C3[6]+C6[3])
    
    C45all = 0.5*(C4[5]+C5[4])
    C46all = 0.5*(C4[6]+C6[4])
    C56all = 0.5*(C5[6]+C6[5])
    
    # Average moduli for cubic crystals
    
    C11cubic = (C11all+C22all+C33all)/3.0
    C12cubic = (C12all+C13all+C23all)/3.0
    C44cubic = (C44all+C55all+C66all)/3.0
    
    bulkmodulus = (C11cubic+2*C12cubic)/3.0
    shearmodulus1 = C44cubic
    shearmodulus2 = (C11cubic-C12cubic)/2.0
    poissonratio = 1.0/(1.0+C11cubic/C12cubic)
    
    printf("========================================= \n")
    printf("Components of the Elastic Constant Tensor \n")
    printf("========================================= \n")
    
    printf("Elastic Constant C11all = %15.8f [GPa] \n", C11all)
    printf("Elastic Constant C22all = %15.8f [GPa] \n", C22all)
    printf("Elastic Constant C33all = %15.8f [GPa] \n", C33all)
    
    printf("Elastic Constant C12all = %15.8f [GPa] \n", C12all)
    printf("Elastic Constant C13all = %15.8f [GPa] \n", C13all)
    printf("Elastic Constant C23all = %15.8f [GPa] \n", C23all)
    
    printf("Elastic Constant C44all = %15.8f [GPa] \n", C44all)
    printf("Elastic Constant C55all = %15.8f [GPa] \n", C55all)
    printf("Elastic Constant C66all = %15.8f [GPa] \n", C66all)
    
    printf("Elastic Constant C14all = %15.8f [GPa] \n", C14all)
    printf("Elastic Constant C15all = %15.8f [GPa] \n", C15all)
    printf("Elastic Constant C16all = %15.8f [GPa] \n", C16all)
    
    printf("Elastic Constant C24all = %15.8f [GPa] \n", C24all)
    printf("Elastic Constant C25all = %15.8f [GPa] \n", C25all)
    printf("Elastic Constant C26all = %15.8f [GPa] \n", C26all)
    
    printf("Elastic Constant C34all = %15.8f [GPa] \n", C34all)
    printf("Elastic Constant C35all = %15.8f [GPa] \n", C35all)
    printf("Elastic Constant C36all = %15.8f [GPa] \n", C36all)
    
    printf("Elastic Constant C45all = %15.8f [GPa] \n", C45all)
    printf("Elastic Constant C46all = %15.8f [GPa] \n", C46all)
    printf("Elastic Constant C56all = %15.8f [GPa] \n", C56all)
    
    printf("========================================= \n")
    printf("Average properties for a cubic crystal    \n")
    printf("========================================= \n")
    
    printf("Bulk Modulus    = %15.8f [GPa] \n", bulkmodulus)
    printf("Shear Modulus 1 = %15.8f [GPa] \n", shearmodulus1)
    printf("Shear Modulus 2 = %15.8f [GPa] \n", shearmodulus2)
    printf("Poisson's ratio = %15.8f [GPa] \n", poissonratio)
}
