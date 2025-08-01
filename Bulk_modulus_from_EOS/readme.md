# Usage

1. Input the lattice constant [Angstrom] and total energy [eV] in `lattice_cohesive_energy`. Select the most stable structure and obtain 5 data points:
   - Variation ratios to test: -6%, -3%, 0%, 3%, 6%
2. Open `fit_eos.gpl` in a text editor and modify the parameters (a, b, c, d). You can use trial and error to find the best fitting values.
3. Run `fit_eos.gpl` (double-click to execute on Windows).
4. The value of “a” found in `fit.log` corresponds to the bulk modulus [eV/Angstrom^3].
   - Tip: 1 [eV/Å³] = 160.2 [GPa]
   - The bulk modulus B is useful for calculating sound velocity using the following formula:

$$
v = \sqrt{ \frac{B}{\rho} }
$$

where $$\( \rho \)$$ is the density of the crystal. Understanding the bulk modulus helps in analyzing mechanical properties of materials.
