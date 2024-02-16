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

```@repl
using Incertus
crd =  CRD()     # complete randomization, targeting 1:1 allocation

w = [1, 2, 3, 4]
crd =  CRD(w)    # complete randomization, targeting allocation specified by w
```

## Permuted Block Design (PBD)

Treatment assignments are made in blocks of size ``bs`` (for a _two-arm_ trial and `1:1` _target_ allocation, ``bs= 2\lambda``; otherwise, ``bs=\lambda W``, where ``W=w_1 + \ldots + w_K`` is a sum of elements of vector ``\mathbf{w}``, _target_ allocation vector. Here, ``\lambda`` is a parameter of the _PBD_, representing the _number of minimal balanced sets in the block of size_ ``bs``).

At the ``j^\text{th}`` allocation step, let ``k^{(j-1)} = \left \lfloor \frac{j-1}{bs}\right\rfloor`` (``\lfloor x \rfloor`` is a `floor` function that returns the greatest integer less than or equal to ``x``). In essence, ``k^{(j-1)}`` is the number of complete blocks among the first ``j-1`` assignments.

The probabilities of treatment assignments within each block are changed according the current imbalance in a block:

- for a _two-arm_ trial and `1:1` _target_ allocation,  
```math
\phi_j = \frac{0.5bs(1+k^{(j-1)})-N_1(j-1)}{bs(1+k^{(j-1)})-(j-1)}, \: j = 1, \ldots, n;
```

- for a _two-arm_ trial with _unequal_ allocation or _multi-arm_ trial with _equal_/_unequal_ allocation,

```math
P_k(j) = \frac{w_k\lambda(1+k^{(j-1)})-N_k(j-1)}{bs(1+k^{(j-1)})-(j-1)}, k = 1, \ldots, K; \: j = 1, \ldots, n.
```
See **[Zhao and Weng (2011), page 955, equation (5)]**

```@docs
PBD
```

```@repl
using Incertus
pbd = PBD(1)      # PBD, targeting 1:1 allocation, with a block size equal to 2*1 = 2
pbd = PBD(3)      # PBD, targeting 1:1 allocation, with a block size equal to 2*3 = 6

w = [1, 2, 3, 4]
pbd =  PBD(w, 1)  # PBD, targeting allocation specified by w, with a block size sum(w)
pbd =  PBD(w, 3)  # PBD, targeting allocation specified by w, with a block size 3*sum(w)
```


## Random Allocation Rule (Rand)

A version of PBD, when the block size ``bs`` equals to the total sample size ``n``. At the ``j^\text{th}`` allocation step, probabilities of treatment assignments are calculated as:

- for a _two-arm_ trial and `1:1` _target_ allocation,  
```math
\phi_j = \frac{0.5n-N_1(j-1)}{n-(j-1)}, \: j = 1, \ldots, n;
```

- for a _two-arm_ trial with _unequal_ allocation or _multi-arm_ trial with _equal_/_unequal_ allocation,

```math
P_k(j) = \frac{nw_k/W-N_k(j-1)}{n-(j-1)}, k = 1, \ldots, K; \: j = 1, \ldots, n,
```
where ``W=w_1 + \ldots + w_K`` is a sum of elements of vector ``\mathbf{w}``, _target_ allocation vector.

```@docs
RAND
```

```@repl
using Incertus

rnd = RAND(50)     # RAND, targeting 1:1 allocation, in a trial with 50 subjects

w = [1, 2, 3, 4]
rnd =  RAND(w, 50) # RAND, targeting allocation specified by w, in a trial with 50 subjects
```


## Truncated Binomial Design (TBD)

Treatment assignments are made with probability 0.5 until one of the treatments receives its quota of ``\frac{n}{2}`` subjects; thereafter all remaining assignments are made deterministically to the opposite treatment.

At the ``j^\text{th}``allocation step, let ``N_1`` and ``N_2`` be the numbers of subjects allocated to treatments s.t. ``N_1+N_2 = j-1.`` Then,

```math
\phi_j = \left\{\begin{array}{rl}
0.5, & \max(N_1, N_2) < \frac{n}{2} \\
1, & N_1 < N2 \\
0, & N_1 > N2
\end{array}\right. ,\: j = 1, \ldots, n.
```

```@docs
TBD
```

```@repl
using Incertus
tbd = TBD(50)  # TBD, targeting 1:1 allocation, in a trial with 50 subjects
```

## Efron's Biased Coin Design (EBCD)

At any allocation step, if treatment numbers ``N_1`` and ``N_2`` are balanced, the next assignment is made with probability 0.5; otherwise, the underrepresented treatment is assigned with probability ``p``, where ``0.5 \leq p \leq 1`` is a fixed and pre-specified parameter that determines the trade-off between balance and randomness.

At the ``j^\text{th}`` allocation step, given treatment numbers ``N_1`` and ``N_2``, s.t. ``N_1+N_2 = j-1``, and imbalance ``d = N_1-N_2``,

```math
\phi_j = \left\{\begin{array}{rl}
0.5, & N_1 = N_2\\
p, & N_1 < N2 \\
1-p, & N_1 > N2
\end{array}\right. ,\: j = 1, \ldots, n.
```

Note that ``p=1`` corresponds to PBD with block size ``b=2``.

```@docs
EBCD
```

```@repl
using Incertus
ebcd = EBCD(2//3) # EBCD, targeting 1:1 allocation, with parameter p=2/3
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

```@repl
using Incertus
abcd = ABCD(2) # ABCD, targeting 1:1 allocation, with parameter a=2
```

## Generalized Biased Coin Design (GBCD)

A generalization of Efron’s BCD. At the ``j^\text{th}`` allocation step, given treatment numbers ``N_1`` and ``N_2``, s.t. ``N_1+N_2 = j-1``, and imbalance ``d = N_1-N_2``,

```math
\phi_j = \left\{\begin{array}{rl}
0.5, & j = 1 \\
\frac{N_2^\gamma}{N_1^\gamma+N_2^\gamma}, & j = 1, \ldots, n.
\end{array}\right. 
```

```@docs
GBCD
```

```@repl
using Incertus
gbcd = GBCD(2) # GBCD, targeting 1:1 allocation, with parameter γ=2
```


## Big Stick Design (BSD)

An example of maximum tolerated imbalance (MTI) procedures. It makes prediction of the future treatment allocations more difficult (even knowing the current sizes of the treatment groups) and controls treatment imbalance at a predefined threshold throughout the experiment. A general MTI procedure specifies a certain boundary for treatment imbalance, say ``mti``, that cannot be exceeded.

At the ``j^\text{th}`` allocation step, given treatment numbers ``N_1`` and ``N_2``, s.t. ``N_1+N_2 = j-1``, and imbalance ``d = N_1-N_2``,

```math
\phi_j = \left\{\begin{array}{rl}
0.5, & |d| < mti \\
0, & d = mti \\
1, & d = -mti 
\end{array}\right. ,\: j = 1, \ldots, n.
```

```@docs
BSD
```

```@repl
using Incertus
bsd = BSD(3) # BSD, targeting 1:1 allocation, with parameter mti=3
```


## Biased Coin Design With Imbalance Tolerance (BCDWIT)

A combination of Efron’s BCD and BSD. At the ``j^\text{th}`` allocation step, given treatment numbers ``N_1`` and ``N_2``, s.t. ``N_1+N_2 = j-1``, and imbalance ``d = N_1-N_2``,

```math
\phi_j = \left\{\begin{array}{rl}
0.5, & |d| < mti\: \& \: d = 0 \\
p, & |d| < mti \: \& \: d < 0 \\
1-p, & |d| < mti \: \& \: d > 0 \\
0, & d = mti \\
1, & d = -mti 
\end{array}\right. ,\: j = 1, \ldots, n.
```

```@docs
BCDWIT
```

```@repl
using Incertus
bcdwit = BCDWIT(2//3, 3) # BCDWIT with p = 2/3 and mti=3, targeting 1:1 allocation
```

## Block Urn Design (BUD)

This design was proposed by **Zhao and Weng (2011)**, to provide a more random design than the PBD.  Let ``N_k(j-1)`` denote the number of treatment ``k`` assignments among first ``j-1`` subjects, and ``k^{(j-1)}=\min\limits_{1 \leq k \leq K} \left \lfloor \frac{N_k(j-1)}{w_k}\right \rfloor`` denote the number of minimal balanced sets among the first ``j-1`` assignments. Then, at the ``j^\text{th}`` allocation step, probabilities of treatment assignments are calculated as:

- for a _two-arm_ trial and `1:1` _target_ allocation,  
```math
\phi_j = \frac{\lambda+\min(N_1(j-1), N_2(j-1))-N_1(j-1)}{2(\lambda+\min(N_1(j-1), N_2(j-1)))-(j-1)}, \: j = 1, \ldots, n;
```

- for a _two-arm_ trial with _unequal_ allocation or _multi-arm_ trial with _equal_/_unequal_ allocation,

```math
P_k(j) = \frac{w_k(\lambda+k^{(j-1)})-N_k(j-1)}{W(\lambda+k^{(j-1)})-(j-1)}, k = 1, \ldots, K; \: j = 1, \ldots, n,
```

where ``W=w_1 + \ldots + w_K`` is a sum of elements of vector ``\mathbf{w}``, _target_ allocation vector.

See **[Zhao and Weng (2011), page 955, equations (2) and (3)]**.

```@docs
BUD
```

```@repl
using Incertus

bud = BUD(2)      # BUD, targeting 1:1 allocation (λ=2)

w = [1, 2, 3, 4]
bud =  BUD(w, 2)  # BUD, targeting allocation specified by w (λ=2)
```


## Ehrenfest Urn Design (EUD)

Another example of the maximum tolerated imbalance (MTI) procedure. At the ``j^\text{th}`` allocation step, given treatment numbers ``N_1`` and ``N_2``, s.t. ``N_1+N_2 = j-1``, and imbalance ``d = N_1-N_2``,

```math
\phi_j = \frac{1}{2}\left(1-\frac{d}{mti}\right) ,\: j = 1, \ldots, n,
```

where ``mti`` (``>0``) is a parameter of the procedure.

```@docs
EUD
```

```@repl
using Incertus
eud = EUD(2) # EUD, targeting 1:1 allocation (mti=2)
```

## Bayesian Biased Coin Design (BBCD)

A special class of _biased coind designs_ (BCDs). At the ``j^\text{th}`` allocation step, given treatment numbers ``N_1`` and ``N_2``, s.t. ``N_1+N_2 = j-1``,

```math
\phi_j = \left\{\begin{array}{rl}
0.5, & j = 1 \\
1, & j = 2 \: \& \: N_1 = 0 \\
0, & j = 2 \: \& \: N_1 = 1 \\
\frac{\left(1 + \frac{N_2}{nN_1}\right)^\frac{1}{\gamma}}{\left(1 + \frac{N_2}{nN_1}\right)^\frac{1}{\gamma} + \left(1 + \frac{N_1}{nN_2}\right)^\frac{1}{\gamma}}, & j \geq 3 
\end{array}\right. ,\: j = 1, \ldots, n.
```

where ``\gamma`` (``>0``) is a parameter of the procedure, and ``n`` is a sample size (pre-specified).

```@docs
BBCD
```

```@repl
using Incertus
bbcd = BBCD(0.1, 40) # BBCD, targeting 1:1 allocation (γ=2, n = 40)
```

# Funcions implemented to calculate allocation probabilities

```@docs
allocation_prb(::CRD)
```

```@docs
allocation_prb(::PBD, ::Vector{Int64})
```

```@docs
allocation_prb(::RAND, ::Vector{Int64})
```

```@docs
allocation_prb(::TBD, ::Vector{Int64})
```

```@docs
allocation_prb(::EBCD, ::Vector{Int64})
```

```@docs
allocation_prb(::ABCD, ::Vector{Int64})
```

```@docs
allocation_prb(::GBCD, ::Vector{Int64})
```

```@docs
allocation_prb(::BSD, ::Vector{Int64})
```

```@docs
allocation_prb(::BCDWIT, ::Vector{Int64})
```

```@docs
allocation_prb(::BUD, ::Vector{Int64})
```

```@docs
allocation_prb(::EUD, ::Vector{Int64})
```

```@docs
allocation_prb(::BBCD, ::Vector{Int64})
```


# Simulation

To perform simulations, the following functionality has been implemented:

```@docs
SimulatedRandomization
```

```@docs
simulate(rnd::Randomization, nsbj::Int64, nsim::Int64, seed::Int64 = 314159)
```

```@docs
simulate(rnd::Vector{Randomization}, nsbj::Int64, nsim::Int64, seed::Int64 = 314159)
```

# Operational characteristics (_two-arm trial with equal 1:1 allocation_)

Several measures of imbalance and randomness have been implemented in the package:

_**Measures of imbalance**_:

- ``\mathbf{E}\left[|D(j)|\right]`` -- expected absolute imbalance at the ``j^\text{th}`` allocation step.
- ``\mathbf{E}\left[|D(j)|^2\right]=\mathbf{var}\left[D(j)\right]`` -- variance of imbalance at the ``j^\text{th}`` allocation step.
- ``\mathbf{E}\left[\max\limits_{1\leq m \leq j}|D(m)|\right]`` -- expected maximum imbalance over the first ``j`` allocation steps.
- ``Imb(j) = \frac{1}{j}\sum\limits_{m=1}^j\frac{\mathbf{E}\left[|D(m)|^2\right]}{m}`` -- a cumulative average loss at the ``j^\text{th}`` allocation step.

_**Measures of randomness**_:

- ``EPCG_{conv}(j) = \frac{1}{j}\sum\limits_{m=1}^j\mathbf{E}\left[G_m = \delta_m\right]`` -- cumulative average of expected proportions of correct guesses over first (``j``) allocation steps under the _convergence_ guessing strategy, where ``G_m`` is a random variable taking values based on the investigator's guess (if ``G_m = \delta_m``, the guess is correct) at the ``m^\text{th}`` allocation step, given current imbalance, ``D(m-1)``:

```math
G_m = \left\{
\begin{array}{rl}
1, & D(m-1) < 0; \\
\sim Bernoulli(0.5), & D(m-1) = 0; \\
0, & D(m-1) > 0.
\end{array}  
\right.  
``` 

- ``EPCG_{max}(j) = \frac{1}{j}\sum\limits_{m=1}^j\mathbf{E}\left[\widetilde{G}_m = \delta_m\right]`` -- cumulative average of expected proportions of correct guesses over first (``j``) allocation steps under the _maximum probability_ guessing strategy, where ``\widetilde{G}_m`` is a random variable taking values based on the investigator's guess (if ``G_m = \delta_m``, the guess is correct) at the ``m^\text{th}`` allocation step, given allocation probability, ``\phi_{m}``:

```math
\widetilde{G}_m = \left\{
\begin{array}{rl}
1, & \phi_m > 0.5; \\
\sim Bernoulli(0.5), & \phi_m = 0.5; \\
0, & \phi_m < 0.5.
\end{array}  
\right.  
```

- ``PD(j) = \frac{1}{j}\sum\limits_{m=1}^j\Pr(\phi_m\in\{0, 1\})`` -- cumulative average of expected proportions of deterministic assignments over first (``j``) allocation steps.
- ``FI(j) = \frac{4}{j}\sum\limits_{m=1}^j\mathbf{E}\left[|\phi_m-0.5|\right]`` -- _forcing index_, which takes values on a scale 0–1. ``FI(j) = 0, \forall j`` for CRD and ``FI(j) = 1`` for PBD with a block size ``bs=2``, assuming $j$ is even (most balanced design). 

_**Balance-randomness trade-off**_:

- ``G(j) = \sqrt{\left\{Imb(j)\right\}^2 + \left\{FI(j)\right\}^2}`` represents a _balance-randomness trade-off_ at the ``j^\text{th}`` allocation step. Lower values of ``G(j)`` indicate better balance–randomness trade-off.

Below, there are functions available for calculating operational characteristics.

## Final imbalance

```@docs
calc_final_imb(sr::SimulatedRandomization)
```

```@docs
calc_final_imb(sr::Vector{SimulatedRandomization})
```

## Expected absolute imbalance

```@docs
calc_expected_abs_imb(sr::SimulatedRandomization)
```

```@docs
calc_expected_abs_imb(sr::Vector{SimulatedRandomization})
```


## Variance of imbalance

```@docs
calc_variance_of_imb(sr::SimulatedRandomization)
```

## Expected maximum imbalance over first allocation steps

```@docs
calc_expected_max_abs_imb(sr::SimulatedRandomization)
```

## A cumulative average loss at the first allocation steps

```@docs
calc_cummean_loss(sr::SimulatedRandomization)
```

## Cumulative average of expected proportions of correct guesses over first allocation steps under _two_ different guessing strategies 

```@docs
calc_cummean_epcg(sr::SimulatedRandomization, gs::String)
```

## Cumulative average of expected proportions of deterministic assignments over first allocation steps

```@docs
calc_cummean_pda(sr::SimulatedRandomization)
```

## Forcing index

```@docs
calc_fi(sr::SimulatedRandomization)
```

## Balance-randomness trade-off

```@docs
calc_brt(sr::SimulatedRandomization)
```


# Visualizing operational characteristics

```@docs
plot(op::DataFrame; kwargs...)
```

```@docs
heatmap(brt::DataFrame; kwargs...)
```

```@docs
violin(final_imb::DataFrame; kwargs...)
```

# Auxiliary functions

```@docs
set_label(::Union{CompleteRandomization, RestrictedRandomization})
```