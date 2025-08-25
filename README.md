Generalized Quadratic Model for Charge Transfer

This library is based on the work published in [Phys. Chem. Chem. Phys](https://doi.org/10.1039/d5cp00866b) **27**, 11318 (2025).

---
# Table of Contents
1. [Reading Outputs](#1-reading-outputs)
    1. [NWChem](#11-nwchem)
2. [Quadratic Interpolations](#2-quadratic-interpolations)
    1. [One-Parabola Model](#21-one-parabola-model)
    2. [Two-Parabola Model](#22-two-parabola-model)
    3. [Generalized Quadratic Model](#23-generalized-quadratic-model)
3. [Charge Transfer](#3-charge-transfer)
    1. [One-Parabola Model](#31-one-parabola-model)
    2. [Two-Parabola Model](#32-two-parabola-model)
    3. [Generalized Quadratic Model](#33-generalized-quadratic-model)
---
# 1. Reading Outputs
Currently we only include a parser for NWChem outputs. But please reach out to us if you are interested in other electronic structure codes, we are happy to help you.

# 1.1 NWChem
To read NWChem outputs, we can use the `read` module from the package. We need to provide the path to the output files for the chemical species.

```python
from cdft import read

cation  = read.nwchem( "path/to/cation/output.out" )
neutral = read.nwchem( "path/to/neutral/output.out" )
anion   = read.nwchem( "path/to/anion/output.out" )
```

Please note that we can use whatever extension we wish for our outputs for as long as we include the complete file name. The `read` module will go through the information in the file and extract the total energy and orbital energies. These attributes are readily accesible from the returned objects.

```python
# Getting the attributes from the 'neutral' parsed output

total_energy = neutral.energy
homo_energy  = neutral.homo
lumo_energy  = neutral.lumo
```

---

# 2. Quadratic Interpolations
We include three different models, namely, the [One-Parabola Model](#21-one-parabola-model) by Parr and Pearson, the [Two-Parabola Model](#22-two-parabola-model) by Gazquez, Cedillo and Vela, and the [Generalized Quadratic Model](#23-generalized-quadratic-model) by Albavera-Mata, Gazquez and Vela.

## 2.1 One-Parabola Model
The One-Parabola Model is based on the assumption that the energy of a system as a function of the number of electrons can be approximated by a single parabola based on the chemical potential $\mu$ and chemical hardness $\eta$. This model is particularly useful for describing systems where the addition or removal of electrons leads to significant changes in energy.

```math
\Delta E(\Delta N) = \mu(\Delta N) + \frac{1}{2} \eta (\Delta N)^2
```

where $`\mu = -(I + A)/2`$ and $`\eta = I - A`$ are defined in terms of the ionization potential $I$ and electron affinity $A$.

```python
from cdft import model

mu, eta = model.one_parabola( cation=cation, neutral=neutral, anion=anion )
```

## 2.2 Two-Parabola Model
The Two-Parabola Model considers two separate parabolas to describe the energy changes associated with electron addition ($\mu^+$ and $\eta^+$) and removal ($\mu^-$ and $\eta^-$). This model provides a more accurate representation for systems where the energy landscape is more complex.

```math
\Delta E^\pm(\Delta N) = \mu^\pm(\Delta N) + \frac{1}{2} \eta (\Delta N)^2
```

where $`\mu^- = -(3I + A)/4`$, $`\mu^+ = -(I + 3A)/4`$, and $`\eta = (I - A)/2`$.

```python
from cdft import model

mu_minus, mu_plus, eta = model.two_parabola( cation=cation, neutral=neutral, anion=anion )
```

## 2.3 Generalized Quadratic Model
The Generalized Quadratic Model further refines the previous models by incorporating information from the highest-occupied and lowest-unnocupied frontier molecular orbitals, HOMO and LUMO, respectively, to better capture the nuances of charge transfer processes. This model is versatile and can be adapted to a wide range of chemical systems.

```math
\Delta E^\pm(\Delta N) = \mu^\pm(\Delta N) + \frac{1}{2} \eta^\pm (\Delta N)^2
```

where $`\mu^- = \varepsilon_\mathrm{HOMO}`$, $`\mu^+ = \varepsilon_\mathrm{LUMO}`$, $`\eta^- = 2(I + \varepsilon_\mathrm{HOMO})`$, and $`\eta^+ = -2(A + \varepsilon_\mathrm{LUMO})`$

```python
from cdft import model

mu_minus, mu_plus, eta_minus_eta_plus = model.generalized( cation=cation, neutral=neutral, anion=anion )
```

---

# 3 Charge Transfer
For the sake of consistency and comparison, we include the same three [quadratic interpolations](#2-quadratic-interpolations). These commonly are used to discern charge trasnfer trends for donor-acceptor reactions, granted that a reactant **A** is donating charge while reactant **B** accepts it. For illustrative purposes, here we will assume that species **B** is our target reagent, where $`{\Delta N}_{A}^{-}`$ and $`{\Delta N}_{A}^{+}`$ denote the nucleophilic and electrophilic charge transfer, respectively.

## 3.1 One-Parabola Model
The charge transfer channels for the One-Parabola Model are given by

```math
\begin{aligned}
    {\Delta N}_{A}^{-} & = +\frac{1}{2} \frac{A_\mathbf{A} - I_\mathbf{B}}{\eta_\mathbf{A} + \eta_\mathbf{B}} \\
    {\Delta N}_{A}^{+} & = -\frac{1}{2} \frac{A_\mathbf{B} - I_\mathbf{A}}{\eta_\mathbf{A} + \eta_\mathbf{B}}
\end{aligned}
```

where $I$, $A$ and $\eta$ are the ionization potential, electron affinity and chemical hardness, respectively, for both interacting species.

```python
from cdft import charge

delta_N_minus, delta_N_plus = charge.one_parabola( cation=cation_A, neutral=neutral_A, anion=anion_A,
                                                   ref_cation=cation_B, ref_neutral=neutral_B, ref_anion=anion_B )
```

## 3.2 Two-Parabola Model
The charge transfer channels for the Two-Parabola Model are given by

```math
\begin{aligned}
    {\Delta N}_{A}^{-} & = \frac{{\mu}^{+}_\mathbf{B} - {\mu}^{-}_\mathbf{A}}{\eta_\mathbf{A} + \eta_\mathbf{B}} \\
    {\Delta N}_{A}^{+} & = \frac{{\mu}^{-}_\mathbf{B} - {\mu}^{+}_\mathbf{A}}{\eta_\mathbf{A} + \eta_\mathbf{B}}
\end{aligned}
```
where $`\mu^\pm`$ and $`\eta`$ are the chemical potential and chemical hardness, respectively, for both interacting species.

```python
from cdft import charge

delta_N_minus, delta_N_plus = charge.one_parabola( cation=cation_A, neutral=neutral_A, anion=anion_A,
                                                   ref_cation=cation_B, ref_neutral=neutral_B, ref_anion=anion_B )
```

## 3.3 Generalized Quadratic Model
The charge transfer channels for the Generalized Quadratic Model are given by

```math
\begin{aligned}
    {\Delta N}_{A}^{-} & = \frac{{\mu}^{+}_\mathbf{B} - {\mu}^{-}_\mathbf{A}}{{\eta}^{-}_\mathbf{A} + {\eta}^{+}_\mathbf{B}} \\
    {\Delta N}_{A}^{+} & = \frac{{\mu}^{-}_\mathbf{B} - {\mu}^{+}_\mathbf{A}}{{\eta}^{+}_\mathbf{A} + {\eta}^{-}_\mathbf{B}}
\end{aligned}
```

where $`\mu^\pm`$ and $`\eta^\pm`$ are the chemical potential and chemical hardness, respectively, for both interacting species.

```python
from cdft import charge

delta_N_minus, delta_N_plus = charge.generalized( cation=cation_A, neutral=neutral_A, anion=anion_A,
                                                  ref_cation=cation_B, ref_neutral=neutral_B, ref_anion=anion_B )
```
