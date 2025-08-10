#!/bin/bash

# Base input file
base_input="case.scf.in"

# Structure optimization using vc-relax
# Output optimized CELL_PARAMETERS and ATOMIC_POSITIONS

# Set number of threads and CPUs
export OMP_NUM_THREADS=1
NCPUs=$(($(nproc) / 2))

# Create log directory if it doesn't exist
mkdir -p log

# Input and output files
input_file="log/case.vcrelax.in"
output_file="log/case.opt.out"

cp "$base_input" "$input_file"

sed -i 's/relax/vc-relax/' "$input_file"

# Run Quantum ESPRESSO vc-relax calculation
mpirun -np ${NCPUs} pw.x < "$input_file" | tee "$output_file"

# Extract optimized CELL_PARAMETERS and ATOMIC_POSITIONS
opt_cell="log/optimized_cell.txt"
opt_atoms="log/optimized_atoms.txt"

# Extract only the numeric rows of CELL_PARAMETERS (3 rows)
grep -A 3 "CELL_PARAMETERS" "$output_file" > "$opt_cell"

# Extract only the coordinates of ATOMIC_POSITIONS (up to the end final coordinates)
awk '/ATOMIC_POSITIONS/,/^[^A-Z]/ {if ($0 ~ /^[A-Z]/) print $0}' "$output_file"  > "$opt_atoms"

echo "Optimized structure saved to $opt_cell and $opt_atoms"

{
  grep "number of atoms/cell" "$output_file" | tail -n 1
  grep "new unit-cell volume" "$output_file" | tail -n 1
  grep "g/cm^3" "$output_file" | tail -n 1
} | tee optimized_data.txt


new_input="case.opt.in"

# Get the value of alat (Bohr)
alat_bohr=$(grep "^CELL_PARAMETERS" "log/optimized_cell.txt" | sed -n 's/.*alat[= ]*\([^)]*\)).*/\1/p')

# Bohr -> Angstrom
alat_angstrom=$(echo "$alat_bohr * 0.529177" | bc -l)

awk -v alat_bohr="$alat_bohr" -v alat_ang="$alat_angstrom" -v cell_file="$opt_cell" -v atom_file="$opt_atoms" '
BEGIN {in_cell=0; in_atoms=0; in_system=0}
{
  if ($0 ~ /^&SYSTEM/) { in_system=1 }
  
  if (in_system && $0 ~ /^\s*A\s*=/) {
    sub(/A\s*=.*/, "A = " alat_ang ",")
  }
  
  if (in_system && $0 ~ /^\/$/) { in_system=0 }
  
  if ($0 ~ /^CELL_PARAMETERS/) {
    while ((getline line < cell_file) > 0) print line
    in_cell=1
    next
  }
  
  if ($0 ~ /^ATOMIC_POSITIONS/) {
    n = 0
    while ((getline line < atom_file) > 0) {
      lines[n++] = line
    }
    for (i = 0; i < n - 1; i++) {
      print lines[i]
    }
    in_atoms = 1
    next
  }
  
  if (in_cell && $0 ~ /^[^ ]/) { in_cell=0 }
  if (in_atoms && $0 ~ /^[^ ]/) { in_atoms=0 }
  
  if (!in_cell && !in_atoms) print $0
}
' "$base_input" > "$new_input"

sed -i '' 's/CELL_PARAMETERS *(alat=[^)]*)/CELL_PARAMETERS {alat}/' case.opt.in

echo "Replaced input saved to $new_input"