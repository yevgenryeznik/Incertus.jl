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