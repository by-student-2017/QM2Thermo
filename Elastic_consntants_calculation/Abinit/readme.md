# Abinit v9.10.4 Installation (Ubuntu 22.04 LTS)
1. sudo apt -y install abinit
2. which abinit

# Shear modulus calculation
1. bash ./run_elastic_ab.sh
2. awk -f compute_elastic_constants_from_stress_ab.awk elastic_results.txt
