# Incertus.jl

The package implements several _restricted randomization_ procedures for _two-_ and _multi-arm_ clinical trials, targeting _equal_ or _unequal allocation_. 

# Randomization, targeting _equal_  `1:1` allocation

The following settings are assumed:

- There are two _treatment arms_ investigated in a trial:
    + ``E`` is an _experimental_ treatment arm.
    + ``C`` is a _control_ treatment arm.
- ``n`` is a total _sample size_.
- ``N_1(j)`` and ``N_2(j)`` are _treatment numbers_, i.e., _sample sizes_ on treatments after the ``j^\text{th}`` allocation step (``N_1(j)+N_2(j) = j, \: j = 1, \ldots, n``).
- ``\delta_j`` is a _treatment indicator_:

```math
\delta_j = \left\{\begin{array}{rl}
1, &\text{if treatment }E \\
0, &\text{if treatment }C
\end{array}\right. .
```

Under these assumptions, a _restricted randomization_ procedure is defined as

```math
\begin{aligned}
\phi_1 &= \Pr(\delta_1 = 1) = 0.5; \\
\phi_j &= \Pr(\delta_j = 1|\delta_1, \ldots, \delta_{j-1}), \: j = 2, \ldots, n.
\end{aligned}
```

# Randomization, targeting _unequal_ allocation in _two-arm_ trial and _equal_/_unequal_ allocation in _multi-arm_ trial

The following settings are assumed:

- ``K`` is a number of treatment arms in a trial (``K \geq 2``).
- ``\mathbf{w} = \left(w_1,\ldots, w_K\right)'`` is the _fixed allocation ratio_, ``w_1:\ldots :w_K``, where ``w_k``'s are positive, not necessarily equal numbers (usually, integers) with the greatest common divisor of 1.
- ``\boldsymbol{\rho}  = \left(\rho_1, \ldots, \rho_K\right)'`` is a vector of _target allocation proportions_, where 
```math
\rho_k = \frac{w_k}{\sum_{k=1}^{K}{w_k}}, 0 \leq \rho_k \leq 1, \text{ and }\sum_{k=1}^K{\rho_k} = 1.
```
- ``n`` is a total _sample size_.
- ``\mathbf{N}(j) = \left(N_1(j), \ldots, N_K(j)\right)'`` is a vector of _treatment numbers_, i.e., numbers of subjects _allocated_ to ``K`` treatments after ``j`` allocations (``1 \leq j \leq n``). 
    + Note that, in general, ``N_k(j)``'s are random variables with ``\sum_{k = 1}^K{N_k(j)} = j``.
- ``\mathbf{P}(j) = \left(P_1(j), \ldots, P_K(j)\right)'`` is a vector of _randomization (allocation) probabilities_ for the ``j^\text{th}``subject. 
    + Note that ``0 \leq P_k(j) \leq 1``, and ``\sum_{k = 1}^K{P_k(j)} = 1`` for each ``j = 1, 2, \ldots, n``. 
    + Also, note that in general, ``\mathbf{P}(j)`` depends on ``\mathbf{N}(j-1)`` (in generalization of Efron's BCD) or on ``\frac{\mathbf{N}(j-1)}{j-1}`` (in generalization of Wei's _urn design_).

Under these assumptions, a _restricted randomization_ procedure is defined as

```math
\begin{aligned}
P_k(1) &= \Pr(\delta_1 = k) = \rho_k, \: k = 1, \ldots, K; \\
P_k(j) &= \Pr(\delta_j = k|\delta_1, \ldots, \delta_{j-1}), \: k = 1, \ldots, K; \: j = 2, \ldots, n.
\end{aligned}
```

# Procedures implemented in the package

## Completely Randomized Design (CRD)

For a _two-arm_ trial and `1:1` _target_ allocation, every subject is allocated to treatments with a fixed probability

```math
\phi_j = 0.5, \: j = 1, \ldots, n.
```

In case of _two-arm_ trial with _unequal_ allocation or _multi-arm_ trial with _equal_/_unequal_ allocation, every subject is allocated to treatments with fixed probabilities that are equal to the target allocation proportions:

```math
P_k(j) = \rho_k, \: k = 1, \ldots, K; \: j = 1, \ldots, n.
```

```@docs
CRD
```

`CRD()` command initializes a _complete_ randomization procedure, targeting `1:1` allocation:

```@repl
using Incertus
crd =  CRD()
```

`CRD(w)` command initializes a _complete_ randomization procedure, targeting allocation specified by `w`:
```@repl
using Incertus # hide
w = [1, 2, 3, 4]
crd =  CRD(w)
```

## Efron's Biased Coin Design (EBCD)

At any allocation step, if treatment numbers ``N_1`` and ``N_2`` are balanced, the next assignment is made with probability 0.5; otherwise, the underrepresented treatment is assigned with probability ``p``, where ``0.5 \leq p \leq 1`` is a fixed and pre-specified parameter that determines the trade-off between balance and randomness.

Note that ``p=1`` corresponds to PBD with block size ``b=2``.

```math
\phi_j = \left\{\begin{array}{rl}
0.5, & N_1 == N_2\\
p, & N_1 < N2 \\
1-p, & N_1 > N2
\end{array}\right. ,\: j = 1, \ldots, n.
```

```@docs
EBCD
```

`EBCD(p)` command initializes a randomization procedure with a parameter ``p``, targeting `1:1` allocation in a trial:

```@repl
using Incertus
ebcd = EBCD(2//3) # a procedure with p=2/3
```

## Adjustable Biased Coin Design (ABCD)

An extension of Efron’s BCD. At the ``j^\text{th}`` allocation step, given treatment numbers ``N_1`` and ``N_2``, s.t. ``N_1+N_2 = j-1``, and imbalance ``d = N_1-N_2``,

```math
\phi_j = \left\{\begin{array}{rl}
0.5, & |d| <= 1 \\
\frac{|d|^a}{1+|d|^a}, & d < -1 \\
\frac{1}{1+|d|^a}, & d > 1 
\end{array}\right. ,\: j = 1, \ldots, n.
```

```@docs
ABCD
```

`ABCD(a)` command initializes a randomization procedure with a parameter ``a``, targeting `1:1` allocation in a trial:

```@repl
using Incertus
abcd = ABCD(2) # a procedure with a=2
```

# Funcions implemented to calculate allocation probabilities

```@docs
allocation_prb(::CRD)
```

```@docs
allocation_prb(::ABCD, ::Vector{Int64})
```

```@docs
allocation_prb(::EBCD, ::Vector{Int64})
```

# Auxiliary functions

```@docs
set_label(::Union{CompleteRandomization, RestrictedRandomization})
```