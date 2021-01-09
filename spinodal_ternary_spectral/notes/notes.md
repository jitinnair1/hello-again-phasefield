---
header-includes:
  - \usepackage[ruled,vlined,linesnumbered]{algorithm2e}
---
# Ternary Coarsening using Cahn-Hilliard Model

For a ternary system, given a specified distribution of precipitates, describe the evolution by computing the composition profile for every `istep` betwen `0` and `nstep`

The free energy used is of the form:

\begin{align}
F=N_{v} \int_{V}\left[f\left(c_{A}, c_{B}, c_{C}\right)+\sum_{i=A, B, C} \kappa_{i}\left(\nabla{c_{i}}\right)^{2}\right] d V
\end{align}

where the bulk component is:

\begin{equation}
{f\left(c_{A}, c_{B}, c_{C}\right)=\frac{1}{2} \sum_{i \neq j} \chi_{i j} c_{i} c_{j}+\sum_{i} c_{i} \ln c_{i}}
\end{equation}

START
  1. Initialise conc_A, conc_B, conc_C, free_energy, dfree_dB, dfree_dC
  2. Create FFTW Plans
  3. Create Initial Microstructure
  4. Print Initial Microstructure to file
  5. Specify Boundary Conditions
  6. Start Time Loop
        * Calculate Free Energy
          - Compute free energy using logarithmic formulation
          - Compute partial derivates of free energy
        * Calculate evolution constants
          - kappa_BB, kappa_CC, kappa_BC
          - mobility_BB, mobility_CC, mobility_BC
        * Move to Fourier Space
        * Update $\tilde{c}$ using governing equation
        * Move to Real Space
        * Print the current value of composition every `iprint` iterations of this loop
  7. Free up memory and cleanup
END
