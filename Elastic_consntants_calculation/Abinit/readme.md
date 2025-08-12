

# Abinit v9.10.4 Installation (Ubuntu 22.04 LTS)
1. sudo apt -y install abinit
2. which abinit

# python
1. pip install numpy

# Shear modulus calculation
1. bash ./run_opt.sh
2. bash ./run_elastic_ab.sh
3. awk -f compute_elastic_constants_from_stress_ab.awk elastic_results.txt
4. python3 compliance_python3.py


# If there is not enough memory
- Reduce the number of parallel processes as follows (in this example it is 1, but if you have enough memory you can use 2 or 4):
```
NCPUs=1
```

# Note
- For conventional cells, values close to those of Materials Projects were obtained. However, for primitive cells, large values were obtained for C14, C56, etc.
- It is better to use relax to optimize the atomic positions. For structures other than the reference structure, calculations will be performed with the settings unchanged, so you can specify relax in case.scf.in or case.opt.in.
- Si Data (elastic.txt) in this code. (Conventional cell calculation. see case.opt.in) (calculation time: 1 h)
```
=========================================
Components of the Elastic Constant Tensor
=========================================
Elastic Constant C11all =    152.07708550 [GPa]
Elastic Constant C22all =    152.07201400 [GPa]
Elastic Constant C33all =    152.07580400 [GPa]
Elastic Constant C12all =     56.18417280 [GPa]
Elastic Constant C13all =     56.18426693 [GPa]
Elastic Constant C23all =     56.18432807 [GPa]
Elastic Constant C44all =     74.62072425 [GPa]
Elastic Constant C55all =     74.62273200 [GPa]
Elastic Constant C66all =     74.62200200 [GPa]
Elastic Constant C14all =      0.00006175 [GPa]
Elastic Constant C15all =      0.00000060 [GPa]
Elastic Constant C16all =     -0.00002153 [GPa]
Elastic Constant C24all =      0.00002343 [GPa]
Elastic Constant C25all =      0.00000105 [GPa]
Elastic Constant C26all =     -0.00002153 [GPa]
Elastic Constant C34all =      0.00002343 [GPa]
Elastic Constant C35all =      0.00000060 [GPa]
Elastic Constant C36all =     -0.00005680 [GPa]
Elastic Constant C45all =      0.00000000 [GPa]
Elastic Constant C46all =      0.00000000 [GPa]
Elastic Constant C56all =      0.00000000 [GPa]
=========================================
Average properties for a cubic crystal
=========================================
Bulk Modulus    =     88.14782657 [GPa]
Shear Modulus 1 =     74.62181942 [GPa]
Shear Modulus 2 =     47.94535595 [GPa]
Poisson's ratio =      0.26978040 [GPa]
```

## Test
- OS: Ubuntu 22.04 LTS (WSL2, Windows 11)
- GPU: 12th Gen Intel(R) Core(TM) i7-12700
- Memory: 32 GB

# Appendix: How do I obtain elastic constants using DFPT?
- Required steps:
1. SCF calculation (normal structure)
2. DFPT response calculation (rfstrs = 1)
3. DDB file generation
4. Analysis with anaddb (elaflag = 1)
5. Extraction of elasticity tensors (Cij)