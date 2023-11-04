# Incertus.jl

The package implements several _restricted randomization_ procedures for _two-_ and _multi-arm_ clinical trials, targeting _equal_ or _unequal allocation_. 

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
allocation_prb(::ABCD, ::Vector{Int64})
```

```@docs
allocation_prb(::EBCD, ::Vector{Int64})
```