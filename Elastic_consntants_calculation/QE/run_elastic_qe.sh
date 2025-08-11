#!/bin/bash

# Strain values to apply
strain_values=(-0.001 +0.001)

# Set number of threads and CPUs
export OMP_NUM_THREADS=1
NCPUs=$(($(nproc) / 2))

# Base input file
#base_input="case.scf.in" # old version
base_input="case.opt.in"  # after "bash run_opt.sh"

# Output file for stress results
results_file="elastic_results.txt"
> "$results_file"
echo "#strain     energy[Ry]      volume[Bohr^3]    s_xx[Ry/Bohr^3] s_xy[Ry/Bohr^3] s_xz[Ry/Bohr^3] s_yy[Ry/Bohr^3] s_yz[Ry/Bohr^3] s_zz[Ry/Bohr^3]" > "$results_file"

# Create log directory if it doesn't exist
mkdir -p log


# set strain = 0 data
dir=0
strain="+0.000"
input_file="log/case.scf.dir${dir}.strain${strain}.in"
output_file="log/case.scf.dir${dir}.strain${strain}.out"

cp "$base_input" "$input_file"

sed -i 's/relax/scf/' "$input_file"

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
printf "%+8.4f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n" \
"$strain" "$energy" "$volume" "$xx" "$xy" "$xz" "$yy" "$yz" "$zz" >> "$results_file"

# Loop over directions (1 to 6)
for dir in {1..6}; do
    
    # Loop over strain values
    for strain in "${strain_values[@]}"; do
        
        A=$(awk '/A / {print $3; exit} /A=/ {print $2; exit}' "$base_input")
        echo "lattice parameter A:", $A
        
        read -r strain <<< $(awk -v strain="${strain}" -v A="${A}" -v dir="${dir}" '
        BEGIN {in_cell=0; line=0}
        /^CELL_PARAMETERS/ {in_cell=1; next}
        in_cell && NF==3 {
          line++
          #---------------------------------------
          for (i=1; i<=3; i++) {
            cell[line,i] = $i * A
          }
          #---------------------------------------
          if (line==3) {
            if (dir == 1) { strain = cell[1,1] * strain }  # e_xx
            if (dir == 2) { strain = cell[2,2] * strain }  # e_yy
            if (dir == 3) { strain = cell[3,3] * strain }  # e_zz
            if (dir == 4) { strain = (cell[2,3] + cell[3,2]) * strain }  # e_yz
            if (dir == 5) { strain = (cell[1,3] + cell[3,1]) * strain }  # e_xz
            if (dir == 6) { strain = (cell[1,2] + cell[2,1]) * strain }  # e_xy
            printf("%+8.4f", strain)
            in_cell=0
            next
          }
          #---------------------------------------
          next
        }
        ' "$base_input")
        
        input_file="log/case.scf.dir${dir}.strain${strain}.in"
        output_file="log/case.scf.dir${dir}.strain${strain}.out"
        
        # Generate strained input file
        awk -v strain="${strain}" -v A="${A}" -v dir="${dir}" '
        BEGIN {in_cell=0; line=0}
        /^CELL_PARAMETERS/ {in_cell=1; print; next}
        in_cell && NF==3 {
          line++
          #---------------------------------------
          for (i=1; i<=3; i++) {
            cell[line,i] = $i * A
          }
          #---------------------------------------
          if (line==3) {
            for (i = 1; i <= 3; i++) {
              for (j = 1; j <= 3; j++) {
                strain_tensor[i, j] = (i == j) ? 1.0 : 0.0
              }
            }
            
            if (dir == 1) { strain_tensor[1,1] += strain }  # e_xx
            if (dir == 2) { strain_tensor[2,2] += strain }  # e_yy
            if (dir == 3) { strain_tensor[3,3] += strain }  # e_zz
            if (dir == 4) { strain_tensor[2,3] += strain }  # e_yz
            if (dir == 5) { strain_tensor[1,3] += strain }  # e_xz
            if (dir == 6) { strain_tensor[1,2] += strain }  # e_xy
            
            for (i = 1; i <= 3; i++) {
              for (j = 1; j <= 3; j++) {
                newcell[i, j] = 0.0
                for (k = 1; k <= 3; k++) {
                  newcell[i, j] += cell[i, k] * strain_tensor[k, j]
                }
              }
              printf("%19.15f %19.15f %19.15f\n", newcell[i,1]/A, newcell[i,2]/A, newcell[i,3]/A)
            }
            in_cell=0
            next
          }
          #---------------------------------------
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
        
        printf "%+8.4f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n" \
        "$strain" "$energy" "$volume" "$xx" "$xy" "$xz" "$yy" "$yz" "$zz" >> "$results_file"
    done
    
    echo "Results for dir=${dir} saved to $results_file"
done

# Evaluate Elastic constants using AWK
echo "command: awk -f compute_elastic_constants_from_stress_qe.awk elastic_results.txt | tee elastic.txt"
awk -f compute_elastic_constants_from_stress_qe.awk elastic_results.txt | tee elastic.txt

python3 compliance_python3.py
