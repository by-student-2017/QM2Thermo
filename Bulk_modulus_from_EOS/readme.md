# Usage

1. Input the lattice constant [Angstrom] and total energy [eV] in `lattice_cohesive_energy`. Select the most stable structure and obtain 5 data points (see https://www.materialscloud.org/discover/sssp or https://molmod.ugent.be/deltacodesdft):
   - Variation ratios to test: -6%, -3%, 0%, 3%, 6%
   - For systems with strong magnetism such as Cr, Mn, Fe, and Co, calculations using spin-polarization are required to obtain elastic moduli close to those obtained experimentally.
   - The calculated values differ from the experimental values by about 20%.
2. Open `fit_eos.gpl` in a text editor and modify the parameters (a, b, c, d). You can use trial and error to find the best fitting values.
3. Run `fit_eos.gpl` (double-click to execute on Windows).
4. The value of “a” found in `fit.log` corresponds to the bulk modulus [eV/Angstrom^3].
   - Tip: 1 [eV/Å³] = 160.2 [GPa]
   - The bulk modulus B is useful for calculating sound velocity using the following formula:


## Sound Velocity Approximation

The sound velocity $$\( v \)$$ can be approximated using:

$$
v_l = \sqrt{ \frac{B + \frac{4}{3}G}{\rho} }, \quad
v_t = \sqrt{ \frac{G}{\rho} }
$$

$$
v_a = \left[ \frac{1}{3} \left( \frac{1}{v_l^3} + \frac{2}{v_t^3} \right) \right]^{-1/3}
$$

$$
\gamma = \frac{9 - 12 \left( \frac{v_t}{v_l} \right)^2}{2 + 4 \left( \frac{v_t}{v_l} \right)^2}
$$

where $$\( \rho \)$$ is the density of the crystal. Understanding the bulk modulus helps in analyzing mechanical properties of materials.


## Elastic Moduli from Bulk Modulus and Poisson's Ratio

Once the bulk modulus $$\( B \)$$ is obtained from EOS fitting, other mechanical properties such as shear modulus $$\( G \)$$ and Young's modulus $$\( E \)$$ can be estimated using the Poisson's ratio $$\( \nu \)$$.

### Key Equations

Given:
- $$\( B \)$$: Bulk modulus [GPa]
- $$\( \nu \)$$: Poisson's ratio (dimensionless)

The following relationships apply:

#### 1. **Shear Modulus** $$\( G \)$$

$$
G = \frac{3B(1 - 2\nu)}{2(1 + \nu)}
$$

#### 2. **Young's Modulus ** $$\( E \)$$

$$
E = 2G(1 + \nu)
$$

Alternatively, you can express $$\( E \)$$ directly from $$\( B \)$$ and $$\( nu \)$$:

$$
E = \frac{9B(1 - \nu)}{(1 + \nu)(1 - 2\nu)}
$$

## Typical Values of Poisson's Ratio

Poisson's ratio $$\( \nu \)$$ varies depending on the material type. It describes the ratio of transverse strain to axial strain when a material is stretched or compressed.

### Common Ranges by Material Type

| Material Type       | Typical $$\( \nu \)$$ Range |
|---------------------|-------------------------|
| Metals              | 0.25 – 0.35             |
| Ceramics            | 0.10 – 0.25             |
| Polymers            | 0.30 – 0.50             |
| Crystalline Solids  | ~0.20 – 0.35            |

- A commonly used default value is **0.3**, especially when experimental or ab initio data is unavailable.
- This value provides a reasonable approximation for many metals, semiconductors, and crystalline materials.

### Notes

- Lower Poisson's ratios (e.g., < 0.2) are typical for brittle materials like ceramics.
- Higher values (e.g., > 0.4) are found in soft, ductile materials such as rubber or certain polymers.
- Accurate Poisson's ratio is essential for calculating shear modulus, Young's modulus, and sound velocities.

