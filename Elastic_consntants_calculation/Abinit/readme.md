# Abinit v9.10.4 Installation (Ubuntu 22.04 LTS)
1. sudo apt -y install abinit
2. which abinit

# python
1. pip install numpy

# Shear modulus calculation
1. bash ./run_elastic_ab.sh
2. awk -f compute_elastic_constants_from_stress_ab.awk elastic_results.txt
3. python3 compliance_python3.py