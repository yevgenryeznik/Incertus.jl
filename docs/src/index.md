# Incertus.jl

The package implements several _restricted randomization_ procedures for _two-_ and _multi-arm_ clinical trials, targeting _equal_ or _unequal allocation_. 

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