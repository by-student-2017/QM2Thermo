#!/bin/bash

# Strain values to apply
strain_values=(-0.010 +0.010)

# Base input file
base_input="case.scf.in"

# output file
results_file="elastic_results.txt"
> "$results_file"
echo "#strain     energy[Ry]      volume[Bohr^3]    s_xx[Ry/Bohr^3] s_xy[Ry/Bohr^3] s_xz[Ry/Bohr^3] s_yy[Ry/Bohr^3] s_yz[Ry/Bohr^3] s_zz[Ry/Bohr^3] Ly[Angstrom]" > "$results_file"

# Set number of threads and CPUs
export OMP_NUM_THREADS=1
NCPUs=$(($(nproc) / 2))


# Create log directory if it doesn't exist
mkdir -p log


# set strain = 0 data
dir=0
strain="+0.000"
input_file="log/case.scf.dir${dir}.strain${strain}.in"
output_file="log/case.scf.dir${dir}.strain${strain}.out"

cp "$base_input" "$input_file"

mpirun -np ${NCPUs} pw.x < "$input_file" | tee "$output_file"

volume=$(awk '/unit-cell volume/ {print $4}' "$output_file")
energy=$(awk '/! *total energy/ {print $5}' "$output_file")

read -r xx xy xz yy yz zz <<< $(awk '
  /Computing stress/ {
  getline; getline; getline;
  printf "%s %s %s ", $1, $2, $3;
  getline;
  printf "%s %s ", $2, $3;
  getline;
  print $3
}' "$output_file")

A=$(awk '/A / {print $3; exit} /A=/ {print $2; exit}' "$base_input")
echo "lattice parameter A: $A"

read Lx0 Ly0 Lz0 <<< $(awk -v A="$A" '
  BEGIN {in_cell=0; line=0}
  /^CELL_PARAMETERS/ {in_cell=1; next}
  in_cell && NF==3 {
    line++
    for (i=1; i<=3; i++) {
      cell[line,i] = $i * A
    }
    if (line==3) {
    # Length
    lx = sqrt(cell[1,1]^2 + cell[1,2]^2 + cell[1,3]^2)
    ly = sqrt(cell[2,1]^2 + cell[2,2]^2 + cell[2,3]^2)
    lz = sqrt(cell[3,1]^2 + cell[3,2]^2 + cell[3,3]^2)
    
    printf "%19.15f %19.15f %19.15f", lx, ly, lz
    exit
    }
  }
' "$base_input")
echo "Lx0: $Lx0, Ly0: $Ly0, Lz0: $Lz0"

printf "%+8.4f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n" \
"$strain" "$energy" "$volume" "$xx" "$xy" "$xz" "$yy" "$yz" "$zz" "$Lx0" "$Ly0" "$Lz0" >> "$results_file"


# Loop over directions (1 to 6)
for dir in {1..6}; do
    
    for strain in "${strain_values[@]}"; do
        input_file="log/case.scf.dir${dir}.strain${strain}.in"
        output_file="log/case.scf.dir${dir}.strain${strain}.out"
        
        A=$(awk '/A / {print $3; exit} /A=/ {print $2; exit}' "$base_input")
        echo "lattice parameter A:", $A
        
        # Generate strained input file
        awk -v strain="${strain}" -v A="${A}" -v dir="${dir}" '
        BEGIN {in_cell=0; line=0}
        /^CELL_PARAMETERS/ {in_cell=1; print; next}
        in_cell && NF==3 {
            line++
            if (line == 1){
              if (dir == 1) {
                $1 = sprintf("%19.15f", $1 + strain/A)
                $2 = sprintf("%19.15f", $2 + strain/A)
                $3 = sprintf("%19.15f", $3 + strain/A)
              }
              if (dir == 5) {
                $3 = sprintf("%19.15f", $3 + strain/A)
              }
              if (dir == 6) {
                $2 = sprintf("%19.15f", $2 + strain/A)
              }
            }
            if (line == 2){
              if (dir == 2) {
                $2 = sprintf("%19.15f", $2 + strain/A)
                $3 = sprintf("%19.15f", $3 + strain/A)
              }
              if (dir == 4) {
                $3 = sprintf("%19.15f", $3 + strain/A)
              }
            }
            if (line == 3){
              if (dir == 3) {
                $3 = sprintf("%19.15f", $3 + strain/A)
              }
            }
            print
            if (line==3) in_cell=0
            next
        }
        {print}
        ' "$base_input" > "$input_file"
        
        mpirun -np ${NCPUs} pw.x < "$input_file" | tee "$output_file"
        
        volume=$(awk '/unit-cell volume/ {print $4}' "$output_file")
        energy=$(awk '/! *total energy/ {print $5}' "$output_file")
        
        read -r xx xy xz yy yz zz <<< $(awk '
            /Computing stress/ {
                getline; getline; getline;
                printf "%s %s %s ", $1, $2, $3;
                getline;
                printf "%s %s ", $2, $3;
                getline;
                print $3
            }' "$output_file")
            
        printf "%+8.4f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n" \
        "$strain" "$energy" "$volume" "$xx" "$xy" "$xz" "$yy" "$yz" "$zz" "$Lx0" "$Ly0" "$Lz0" >> "$results_file"
    done
    
    echo "Results for dir=${dir} saved to $results_file"
    echo "command: awk -f compute_shear_modulus_from_stress_qe.awk $results_file"
    awk -f compute_shear_modulus_from_stress_qe.awk "$results_file"
done

# Evaluate Elastic constants using AWK
echo "command: awk -f compute_elastic_constants_from_stress_qe.awk elastic_results.txt"
awk -f compute_elastic_constants_from_stress_qe.awk elastic_results.txt

python3 compliance_python3.py
