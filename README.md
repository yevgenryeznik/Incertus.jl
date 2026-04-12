# Incertus.jl

# Incertus.jl

> *"The Julia Lego Blocks for Randomized Clinical Trial Designs"*

A `Julia` package to simulate randomization procedures for **two- and multi-arm clinical trials** targeting **equal or unequal allocation**. Incertus.jl generates treatment randomization sequences of a given length and evaluates the operating characteristics of any chosen procedure through Monte Carlo simulation — fast enough to use interactively, and open-ended enough to accommodate new designs.

[![v0.1.0](https://img.shields.io/badge/release-v0.1.0-blue)](https://github.com/yevgenryeznik/Incertus.jl/releases/tag/v0.1.0)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-green)](LICENSE)
[![Julia](https://img.shields.io/badge/language-Julia-purple)](https://julialang.org/)
[![Docs](https://img.shields.io/badge/docs-online-orange)](https://yevgenryeznik.github.io/Incertus.jl/build/)

---

## Overview

Randomization is the cornerstone of any comparative clinical trial. The choice of randomization procedure involves a fundamental trade-off between two competing objectives:

- **Balance** — keeping treatment group sizes close to the target allocation ratio throughout the trial
- **Randomness** — making the sequence of treatment assignments unpredictable, thereby guarding against selection bias

Incertus.jl provides a unified, simulation-based framework for exploring this trade-off across a broad catalogue of randomization designs, for trials with any number of arms and any target allocation ratio — including irrational-valued proportions that arise from optimal design theory.

### Associated Paper

> Ryeznik, Y. and Sverdlov, O. *Incertus.jl — The Julia Lego Blocks for Randomized Clinical Trial Designs.* arXiv preprint, 2024. [arXiv:2407.14248](https://arxiv.org/abs/2407.14248)

### Documentation

Full package documentation is available at:
**[yevgenryeznik.github.io/Incertus.jl/build/](https://yevgenryeznik.github.io/Incertus.jl/build/)**

---

## Features

- **Multi-arm support** — works for K ≥ 2 treatment arms with any fixed target allocation ratio ρ₁ : ρ₂ : … : ρ_K
- **Equal and unequal allocation** — both integer-ratio (e.g., 1:1, 1:2, 2:3:1) and irrational-valued optimal allocation proportions
- **Large procedure catalogue** — implements the most important restricted randomization designs from the clinical trial literature, covering the full balance/randomness spectrum
- **Monte Carlo simulation engine** — evaluates statistical operating characteristics for any design at any sample size
- **R interoperability** — can be called from R via `JuliaCall` or `JuliaConnectoR`
- **Open-ended architecture** — new randomization procedures can be plugged in alongside the existing ones
- **Validation tool** — useful for verifying randomization methods for which no other software implementation is readily available

---

## Randomization Procedures Implemented

### Group 1 — Designs for Integer-Ratio (Equal or Simple Unequal) Allocation

These procedures are suited for trials where the target ratio can be expressed as small positive integers (e.g., 1:1, 2:1, 1:2:1). Many enforce a **Maximum Tolerated Imbalance (MTI)** constraint — a pre-specified bound on how far the running allocation can deviate from target at any point in the trial.

| Abbreviation | Procedure | Key Property |
|---|---|---|
| **CRD** | Complete Randomization Design | No imbalance control; maximum randomness |
| **PBD** | Permuted Block Design | Perfect balance within blocks; highest predictability |
| **BSD** | Big Stick Design (Soares & Wu, 1983) | MTI; forces balance only when bound is hit |
| **EBCD** | Efron's Biased Coin Design (Efron, 1971) | Probabilistic balance control; tunable via bias parameter *p* |
| **ABCD** | Accelerated Biased Coin Design (Wei, 1977) | Urn-based; balance and randomness trade-off |
| **BUD** | Block Urn Design (Zhao & Weng, 2011) | Urn-based; resets blocks to reduce predictability |
| **DLUD** | Drop-the-Loser Urn Design (Ivanova, 2003) | Low variability urn design; near-optimal balance |

### Group 2 — Designs for Any Fixed Allocation (Including Irrational Proportions)

These procedures can target any fixed allocation proportion ρ₁ : … : ρ_K, including non-integer (irrational) values arising from optimal design theory (e.g., D-optimal allocation in dose-response studies).

| Abbreviation | Procedure | Key Property |
|---|---|---|
| **CRD** | Complete Randomization (multinomial) | Assigns independently with P(arm k) = ρ_k |
| **MWUD** | Mass Weighted Urn Design (Zhao, 2015) | Controlled maximum imbalance; good balance/randomness trade-off |
| **DBCD** | Doubly Adaptive Biased Coin Design (Hu & Zhang, 2004) | Adaptive; targets any fixed allocation with high precision |
| **MaxEnt** | Maximum Entropy Constrained Balance Randomization (Klotz, 1978) | Maximises allocation entropy subject to a balance constraint; tunable via η |

---

## Performance Metrics

All procedures are evaluated through the following operating characteristics, estimated via Monte Carlo simulation:

| Metric | Symbol | What It Measures |
|--------|--------|-----------------|
| **Maximum Probability of Maximum Imbalance** | MPM(n) | Expected worst-case deviation from target allocation across all steps up to n |
| **Average Squared Deviation** | ASD(n) | Average squared imbalance from target allocation; measures overall balance |
| **Forcing Index** | FI(n) | Proportion of allocation steps that are deterministic; measures loss of randomness |
| **Selection Bias Factor** | SBF | Expected proportion of assignments correctly guessed by an unblinded observer |
| **Overall Distance** | D(n) | Combined balance/randomness summary: Euclidean distance from the (0,0) ideal in (ASD, FI) space |
| **ARP Property** | — | Allocation Ratio Preserving: whether unconditional treatment probabilities equal ρ_k at every step |

The ARP property is assessed via a plot of unconditional randomization probabilities π_{jk} against the allocation step j, for each arm k. Departures from ρ_k signal a non-ARP procedure.

---

## Installation

From the Julia package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/yevgenryeznik/Incertus.jl")
```

Or in `Pkg` REPL mode (press `]`):

```
pkg> add https://github.com/yevgenryeznik/Incertus.jl
```

---

## Quick Start

```julia
using Incertus

# --- 1. Define a target allocation for a 3-arm trial (1:1:1) ---
n      = 60
target = [1, 1, 1]   # equal allocation

# --- 2. Generate a single randomization sequence using a Permuted Block Design ---
seq = pbd(n, target, block_size=6)

# --- 3. Evaluate operating characteristics via Monte Carlo ---
nsim    = 10_000
results = simulate(PBD(block_size=6), n, target; nsim=nsim)
# results contains MPM, ASD, FI, SBF, and ARP estimates at each allocation step
```

For unequal allocation with irrational proportions (e.g., D-optimal for a 3-arm dose-response study):

```julia
target_proportions = [0.407, 0.336, 0.257]   # D-optimal allocation proportions

results_mwud = simulate(MWUD(alpha=10), n, target_proportions; nsim=10_000)
results_dbcd = simulate(DBCD(gamma=2),  n, target_proportions; nsim=10_000)
```

---

## Repository Structure

```
Incertus.jl/
├── src/              # Package source: procedure implementations, simulation engine, metrics
├── docs/             # Documentation source (Documenter.jl)
├── Project.toml      # Package dependencies and metadata
├── Manifest.toml     # Resolved dependency versions
├── LICENSE           # GPL-3.0
└── README.md
```

---

## Dependencies

All dependencies are managed via Julia's built-in package manager:

| Package | Purpose |
|---------|---------|
| `Distributions` | Probability distributions for sequence generation |
| `DataFrames` | Tabular output of simulation results |
| `Plots` + `StatsPlots` | Visualization of balance/randomness profiles and ARP plots |
| `ColorSchemes` | Plot aesthetics |
| `Latexify` | LaTeX-formatted output of randomization sequences and summary tables |
| `ProgressMeter` | Progress display for long Monte Carlo runs |
| `Roots` | Root-finding for design parameter calibration |
| `Statistics` + `StatsBase` | Summary statistics for simulation outputs |
| `Match` | Pattern matching for procedure dispatch |
| `Pipe` | Pipe operator for readable data processing chains |

---

## R Interoperability

Incertus.jl can be called from R via the `JuliaCall` package:

```r
library(JuliaCall)
julia_setup()
julia_library("Incertus")

# Generate a PBD sequence from R
julia_eval('pbd(60, [1, 1, 1], block_size=6)')
```

---

## Authors

| Name | Affiliation |
|------|-------------|
| [Yevgen Ryeznik](https://github.com/yevgenryeznik) | Department of Mathematics, Uppsala University |
| Oleksandr Sverdlov | Early Development Analytics, Novartis Pharmaceuticals |

---

## License

Released under the [GNU General Public License v3.0](LICENSE).

---

## Citation

If you use Incertus.jl in your work, please cite:

```bibtex
@misc{ryeznik2024incertus,
  author        = {Ryeznik, Yevgen and Sverdlov, Oleksandr},
  title         = {{Incertus.jl} -- The {Julia} {Lego} Blocks for Randomized Clinical Trial Designs},
  year          = {2024},
  eprint        = {2407.14248},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME},
  url           = {https://arxiv.org/abs/2407.14248}
}
```
