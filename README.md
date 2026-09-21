# 1D Stationary Schrödinger Equation — Scattering States

Numerical solution of the one-dimensional, time-independent Schrödinger equation for a particle hitting a
potential barrier. The program computes the real and imaginary parts of the wavefunction ψ(x), and the
reflection (R) and transmission (T) coefficients, using two independent numerical methods. For a simple
rectangular barrier it checks both against the analytical result.

> M1 Physics computing project — Université Claude Bernard Lyon 1 (April 2023)
> Authors: Corentin Le Carrer, Aras Selahiye — Supervisor: Abdul Rahman Allouche

## Physics

$$-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2} + V(x)\,\psi(x) = E\,\psi(x)$$

Atomic units are used (ħ = 1).

- **Method 1: finite-difference integration.** Start from a pure transmitted plane wave on the right of the
  barrier and integrate the discretised equation backwards through the potential. Then match ψ and ψ′ to
  incident + reflected waves on the left to get R and T.
- **Method 2: minimisation.** Find the wave amplitudes by minimising the boundary-condition mismatch with
  a steepest-descent algorithm (numerical gradient).

The program also scans T as a function of energy, and reports the first resonance energy, its transmission
and its position.

## Available potentials

| # | Potential |
|---|-----------|
| 1 | Rectangular barrier |
| 2 | Gaussian barrier |
| 3 | Double barrier (resonant tunnelling) |
| 4 | Triple barrier |
| 5 | Triangular barrier |
| 6 | Alpha-decay (Coulomb-like) barrier |

## Repository layout

```
Schrödinger ProgramFile/   C++ source (main.cpp, schrodinger.cpp/.h), Makefile, input file valeur.txt
Differents Résultats/      Example outputs (double barrier) and older result plots (PDF)
Rapport/RapportModNum.pdf  Full project report (French)
```

## Build

```bash
cd "Schrödinger ProgramFile"
make            # produces ./sch   (a prebuilt Windows sch.exe is also included)
```

Requires a C++17 compiler (g++).

## Usage

1. Edit `valeur.txt`. Its first lines are read in this order:

   | Line | Meaning |
   |------|---------|
   | 1 | Output directory to create (path + folder name) |
   | 2 | Potential type (1–6, see table above) |
   | 3 | `E` — particle energy |
   | 4 | `m` — particle mass |
   | 5 | `a` — barrier width |
   | 6 | `N` — number of points for ψ |
   | 7 | `V0` — potential height |
   | 8 | `lambda` — extra shape parameter |

2. Run `./sch` from the same directory. The output folder is created at the given path.

## Output files

| File | Content |
|------|---------|
| `ComplexPsiValues1/2/3.txt` | ψ(x) before / inside / after the barrier — columns: `x, Re ψ (m1), Im ψ (m1), Re ψ (m2), Im ψ (m2)` |
| `potentiel.txt` | V(x) over the barrier width |
| `Coeff_transmition_fonction_Eeta.txt` | T(E) and T(x) |
| `resultat.txt` | Input parameters, numerical vs. analytical R and T for both methods, their deviation, first resonance |

## Plotting (gnuplot)

```gnuplot
plot "ComplexPsiValues2.txt" u 1:2 w l   # Re ψ, method 1
plot "ComplexPsiValues2.txt" u 1:3 w l   # Im ψ, method 1
plot "ComplexPsiValues2.txt" u 1:4 w l   # Re ψ, method 2
plot "ComplexPsiValues2.txt" u 1:5 w l   # Im ψ, method 2
```

## License

[Mozilla Public License 2.0](LICENSE)
