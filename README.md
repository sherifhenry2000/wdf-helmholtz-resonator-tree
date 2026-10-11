# Coupled Helmholtz Resonator Tree via Wave Digital Filters (WDF)

A discrete-time physical model of a coupled acoustic Helmholtz resonator tree simulated in MATLAB using **Wave Digital Filters (WDF)**. The network translates continuous-time lumped electro-acoustic models into discrete wave-scattering networks via the bilinear transform, guaranteeing structural passivity and zero delay-free loops without matrix inversion at each time step.

The model is verified against an analog-equivalent physical **Simscape (SSC)** ground truth, achieving a Normalized Mean Squared Error (**NMSE**) of **$2.489 \times 10^{-25}$**.

---

## Electro-Acoustic Lumped Modeling

Using the **impedance analogy**, acoustic pressure $p(t)$ maps to electrical voltage $v(t)$, and volume velocity $U(t)$ maps to current $i(t)$:

$$p(t) \longleftrightarrow v(t), \quad U(t) \longleftrightarrow i(t)$$

Each resonator cavity and interconnecting neck is modeled as lumped circuit elements:
- **Acoustic Inertance / Mass (Neck):** Modeled as an inductor $L = \frac{\rho \cdot l}{S}$
- **Acoustic Compliance (Cavity):** Modeled as a capacitor $C = \frac{V}{\rho \cdot c^2}$
- **Acoustic Losses (Radiation & Viscosity):** Modeled as a resistor $R = \frac{\rho \cdot c}{S}$

### System Parameters ($f_s = 96\text{ kHz}$)
- **Air Density ($\rho$):** $1.2\text{ kg/m}^3$
- **Speed of Sound ($c$):** $343.2\text{ m/s}$
- **Volumes ($V_{00}, V_{10}, V_{11}, V_{20}, V_{21}$):** $0.10\text{ m}^3$
- **Neck Lengths ($l_{00}, l_{10}, l_{11}, l_{20}, l_{21}$):** $10.0\text{ m}$
- **Neck Areas ($S_{00}, S_{10}, S_{11}, S_{20}, S_{21}$):** $100.0\text{ m}^2$

---

## Wave Digital Formulation

Kirchhoff variables map to discrete wave quantities via port resistance $Z[k]$:

$$a[k] = v[k] + Z[k] \cdot i[k], \quad b[k] = v[k] - Z[k] \cdot i[k]$$

### Element Port Resistances
Applying the bilinear transform with sampling period $T_s = 1/f_s$:
- **Resistors (Memoryless Leaves):** $b[k] = 0 \implies Z_R = R$
- **Inductors (Delay Leaves):** $b[k] = -a[k-1] \implies Z_L = \frac{2L}{T_s}$
- **Capacitors (Delay Leaves):** $b[k] = a[k-1] \implies Z_C = \frac{T_s}{2C}$
- **Ideal Voltage Source (Root):** $b[k] = 2V_{\text{in}}[k] - a[k]$ (Unadapted)

### Multi-Port Junction Adaptors
Scattering is solved per adaptor using vectorized scattering formulations:
- **Series Adaptors ($N$-Port):**
  $$\mathbf{S} = \mathbf{I} - \boldsymbol{\alpha}^T \mathbf{1}_{1 \times N}, \quad \alpha_j = \frac{2 Z_j}{\sum_{k=1}^N Z_k}$$
  *Port 1 reflection-free condition:* $Z_1 = \sum_{k=2}^N Z_k$

- **Parallel Adaptors ($N$-Port):**
  $$\mathbf{S} = \mathbf{1}_{N \times 1} \boldsymbol{\alpha} - \mathbf{I}, \quad \alpha_j = \frac{2 G_j}{\sum_{k=1}^N G_k}, \quad G_j = \frac{1}{Z_j}$$
  *Port 1 reflection-free condition:* $Z_1 = \left( \sum_{k=2}^N Z_k^{-1} \right)^{-1}$

### Tree Adaptation Hierarchy
1. **Peripheral Series Adaptors:**
   - Adaptor 25: $Z_{25} = Z_{28} + Z_{26} + Z_{27} \quad (R_{20}, L_{20}, C_{20})$
   - Adaptor 21: $Z_{21} = Z_{22} + Z_{24} + Z_{23} \quad (R_{21}, L_{21}, C_{21})$
   - Adaptor 9: $Z_9 = Z_{10} + Z_{12} + Z_{11} \quad (R_{11}, L_{11}, C_{11})$
2. **Intermediate Cavity Parallel Adaptor:**
   - Adaptor 17: $Z_{18} = Z_{25}, \quad Z_{20} = Z_{21}$
   - $Z_{17} = \left( Z_{18}^{-1} + Z_{19}^{-1} + Z_{20}^{-1} \right)^{-1} \quad (C_{10})$
3. **Mid-Branch Series Adaptor:**
   - Adaptor 13: $Z_{15} = Z_{17}, \quad Z_{13} = Z_{16} + Z_{14} + Z_{15} \quad (R_{10}, L_{10})$
4. **Primary Cavity Parallel Adaptor:**
   - Adaptor 5: $Z_8 = Z_9, \quad Z_7 = Z_{13}$
   - $Z_5 = \left( Z_6^{-1} + Z_7^{-1} + Z_8^{-1} \right)^{-1} \quad (C_{00})$
5. **Root Series Adaptor:**
   - Adaptor 1: $Z_3 = Z_5, \quad Z_1 = Z_2 + Z_4 + Z_3 \quad (R_{00}, L_{00})$


