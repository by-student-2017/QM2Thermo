#!/bin/bash

# Set number of threads and CPUs
export OMP_NUM_THREADS=1
NCPUs=$(($(nproc) / 2))

# Base input file
base_input="case.scf.in"

# Output file for stress results
results_file="shear_results.txt"
> "$results_file"
echo "#strain     energy[eV]      volume[bohr^3]    s_xx[GPa]       s_xy[GPa]       s_xz[GPa]       s_yy[GPa]       s_yz[GPa]       s_zz[GPa] Ly[Angstrom]" > "$results_file"

# Strain values to apply
strain_values=(-0.010 -0.005 +0.000 +0.005 +0.010)

# Create log directory if it doesn't exist
mkdir -p log

# Loop over strain values
for strain in "${strain_values[@]}"; do
    input_file="log/case.scf.${strain}.in"
    output_file="log/case.scf.${strain}.out"
    
    # Get A in &SYSTEM section
    A=$(awk '/acell / {print $2; exit}' "$base_input")
    B=$(awk '/acell / {print $3; exit}' "$base_input")
    C=$(awk '/acell / {print $4; exit}' "$base_input")
    echo "lattice parameter A: ${A}"
    echo "lattice parameter B: ${B}"
    echo "lattice parameter C: ${C}"
    
    # Generate strained input file using awk
    awk -v strain="${strain}" -v A="$A" -v B="$B" -v C="$C" '
    BEGIN {in_cell=0; line=0}
    /rprim/ {in_cell=1}
    in_cell {
        line++
        #---------------------------------------
        # Distortion is introduced in this range.
        if (line==2) {
            $1 = sprintf("%19.15f", $1 + (strain/B))
        }
        #---------------------------------------
        print
        if (line==3) in_cell=0
        next
    }
    {print}
    ' "$base_input" > "$input_file"

    # Get Ly
    Ly=$(awk -v A="$A" -v B="$B" -v C="$C" '
    BEGIN {in_cell=0; line=0}
    /rprim/ {in_cell=1}
    in_cell {
        line++
        #---------------------------------------
        # Distortion is introduced in this range.
        if (line==2) { # a2 line
          len = sqrt(($1*A)*($1*A) + ($2*B)*($2*B) + ($3*C)*($3*C))
          printf("%19.15f", len)
          exit
        }
        #---------------------------------------
        if (line==3) in_cell=0
        next
    }
    ' "$base_input")
    echo "lattice parameter Ly:", ${Ly}

    # Run Abinit and extract stress tensor
    mpirun -np ${NCPUs} abinit "$input_file" | tee "$output_file"

    # Extract unit-cell volume
    volume=$(awk '/Unit cell volume ucvol=/ {print $5}' "$output_file")

    # Extract total energy
    energy=$(awk '/Etot       = :/ {print $7}' "$output_file")

    # Extract all 6 components of the stress tensor (Ry/Bohr^3)
    read -r xx yz yy xz zz xy <<< $(awk '
        /-Cartesian components of stress tensor \(GPa\)/ {
            getline;
            printf "%s %s ", $4, $7;
            getline;
            printf "%s %s ", $4, $7;
            getline;
            printf "%s %s ", $4, $7;
        }' "$output_file")

    # Output strain, energy, volume, and stress tensor components
    printf "%+8.4f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\\n" \
    "$strain" "$energy" "$volume" "$xx" "$xy" "$xz" "$yy" "$yz" "$zz" "$Ly">> "$results_file"

done

echo "Shear strain calculations completed. Results saved to $results_file."
echo ""
# Evaluate shear modulus using AWK
echo "command: awk -f compute_shear_modulus_from_stress_ab.awk shear_results.txt"
awk -f compute_shear_modulus_from_stress_ab.awk shear_results.txt

rm -f caseo_DDB caseo_DEN caseo_EBANDS.agr caseo_WFK
rm -f caseo_EIG caseo_EIG.nc caseo_GSR.nc caseo_OUT.nc 
