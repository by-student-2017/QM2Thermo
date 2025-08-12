#!/bin/bash

# Base input file
base_input="case.scf.in"

# Set number of threads and CPUs
export OMP_NUM_THREADS=1
NCPUs=$(($(nproc) / 2))

mkdir -p log

input_file="log/case.opt.in"
output_file="log/case.opt.out"

cp "$base_input" "$input_file"

sed -i '/^[[:space:]]*ionmov[[:space:]]/s/.*/ionmov 2/' "$input_file"
sed -i '/^[[:space:]]*optcell[[:space:]]/s/.*/optcell 2/' "$input_file"

# Run Abinit cell+atom optimization
mpirun -np ${NCPUs} abinit "$base_input" | tee "$output_file"

# Extract optimized CELL_PARAMETERS and ATOMIC_POSITIONS
opt_cell="log/optimized_cell.txt"
opt_atoms="log/optimized_atoms.txt"

natom=$(grep "natom" "$output_file" | tail -1 | awk '{print $2}')
natom_one_minus=$(grep "natom" "$output_file" | tail -1 | awk '{print $2-1}')
echo "natom: ${natom}"

acell=$(grep "acell" "$output_file" | tail -1 | awk '{$5="# "$5; print}')

# Extract only the numeric rows of  (3 rows)
grep -A 2 "rprim" "$output_file" | tail -3 > "$opt_cell"

# Extract only the coordinates of ATOMIC_POSITIONS (up to the end final coordinates)
grep -A "$natom_one_minus" "xred" "$output_file" | tail -"${natom}" > "$opt_atoms"

new_input="case.opt.in"

awk -v acell="$acell" -v cell_file="$opt_cell" -v atom_file="$opt_atoms" '
BEGIN {in_cell=0; in_atoms=0}
{
  if ($0 ~ /acell /) {
    sub(/acell.*/, acell)
  }
  
  if ($0 ~ /rprim /) {
    while ((getline line < cell_file) > 0) print line
    print ""
    in_cell=1
    next
  }
  
  if ($0 ~ /xred /) {
    n = 0
    while ((getline line < atom_file) > 0) {
      lines[n++] = line
    }
    for (i = 0; i < n; i++) {
      print lines[i]
    }
    print ""
    in_atoms = 1
    next
  }
  
  if (in_cell && $0 ~ /^[^ ]/) { in_cell=0 }
  if (in_atoms && $0 ~ /^[^ ]/) { in_atoms=0 }
  
  if (!in_cell && !in_atoms) print $0
}
' "$base_input" > "$new_input"

echo "Replaced input saved to $new_input"

rm -f case.scf.abo caseo_DDB caseo_EBANDS.agr 
rm -f caseo_EIG caseo_EIG.nc caseo_GSR.nc caseo_HIST.nc
rm -f caseo_OUT.nc caseo_TIM*_DEN caseo_WFK
