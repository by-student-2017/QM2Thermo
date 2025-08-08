#!/bin/bash

# Set number of threads and CPUs
export OMP_NUM_THREADS=1
NCPUs=$(($(nproc) / 2))

# Base input file
base_input="case.scf.in"

# Output file for stress results
results_file="bulk_results.txt"
> "$results_file"
echo "#strain     energy[Ry]      volume[Bohr^3]    s_xx[Ry/Bohr^3] s_xy[Ry/Bohr^3] s_xz[Ry/Bohr^3] s_yy[Ry/Bohr^3] s_yz[Ry/Bohr^3] s_zz[Ry/Bohr^3]" > "$results_file"

# Strain values to apply (Note: strain_values/A*100 [%], A in &SYSTEM section)
#strain_values=(-0.012 -0.011 -0.010 -0.009 -0.008 -0.007 -0.006 -0.005 -0.004 -0.003 -0.002 -0.001
#        +0.000 +0.001 +0.002 +0.003 +0.004 +0.005 +0.006 +0.007 +0.008 +0.009 +0.010 +0.011 +0.012)
#--------------
#strain_values=(-0.012 -0.010 -0.008 -0.006 -0.004 -0.002
#        +0.000 +0.002 +0.004 +0.006 +0.008 +0.010 +0.012)
#--------------
#strain_values=(-0.012 -0.009 -0.006 -0.003
#        +0.000 +0.003 +0.006 +0.009 +0.012)
#--------------
#strain_values=(-0.012 -0.008 -0.004
#        +0.000 +0.004 +0.008 +0.012)
#--------------
strain_values=(-0.020 -0.015 -0.010 -0.005
        +0.000 +0.005 +0.010 +0.015 +0.020)
#--------------

# Create log directory if it doesn't exist
mkdir -p log

# Loop over strain values
for strain in "${strain_values[@]}"; do
    input_file="log/case.scf.${strain}.in"
    output_file="log/case.scf.${strain}.out"

    # Get A in &SYSTEM section
    A=$(awk '/A / {print $3; exit} /A=/ {print $2; exit}' "$base_input")
    #A=$(awk '
    #  /A / {print $3; exit}
    #  /A=/ {print $2; exit}
    #' "$base_input")
    echo "lattice parameter A: ${A}"

    # Generate strained input file using awk
    awk -v strain="${strain}" -v A="${A}" '
    BEGIN {in_cell=0; line=0}
    /^CELL_PARAMETERS/ {in_cell=1; print; next}
    in_cell && NF==3 {
        line++
        #---------------------------------------
        # Distortion is introduced in this range.
        for (i=1; i<=3; i++) {
            $i = sprintf("%19.15f", $i * (1 + strain/A))
        }
        #---------------------------------------
        print
        if (line==3) in_cell=0
        next
    }
    {print}
    ' "$base_input" > "$input_file"
    
    # Run QE and extract stress tensor
    mpirun -np ${NCPUs} pw.x < "$input_file" | tee "$output_file"

    # Extract unit-cell volume
    volume=$(awk '/unit-cell volume/ {print $4}' "$output_file")

    # Extract total energy
    energy=$(awk '/! *total energy/ {print $5}' "$output_file")

    # Extract all 6 components of the stress tensor (Ry/Bohr^3)
    read -r xx xy xz yy yz zz <<< $(awk '
        /Computing stress/ {
            getline; getline; getline;
            printf "%s %s %s ", $1, $2, $3;
            getline;
            printf "%s %s ", $2, $3;
            getline;
            print $3
        }' "$output_file")
    # Output strain, energy, volume, and stress tensor components
    printf "%+8.4f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\\n" \
    "$strain" "$energy" "$volume" "$xx" "$xy" "$xz" "$yy" "$yz" "$zz" >> "$results_file"

done

gnuplot fit_eos_qe.gpl
