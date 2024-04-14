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


## Mass Weighted Urn Design (MWUD)

A randomization procedure that can be used fo for two- and multi-arm trials, targeting both equal and unequal allocation. At the ``j^\text{th}`` allocation step, given ``K`` treatments and corresponding treatment numbers ``N_1(j-1) \ldots N_K(j-1)``, s.t. ``\sum_{k = 1}^K{N_k(j-1)} = j-1``,

```math
P_k(j) = \frac{\max\left\{\alpha\rho_k-N_k(j-1) +(j-1)\rho_k , 0\right\}}{\sum_{k=1}^K\max\left\{\alpha\rho_k-N_k(j-1) +(j-1)\rho_k , 0\right\}} ,\: j = 1, \ldots, n.
```

where ``\alpha`` (``>0``) is a parameter of the procedure that controls _maximum tolerated imbalance_.

```@docs
MWUD
```

```@repl
using Incertus
mwud = MWUD(2)     # MWUD, targeting 1:1 allocation (α=2)

w = [1, 2, 3, 4]   # target allocation ratio
mwud =  MWUD(w, 2) # MWUD, targeting 1:2:3:4 allocation (α=2)
```


## Doubly-Adaptive Biased Coin Design (DBCD)

A randomization procedure that can be used fo for two- and multi-arm trials with ``K`` treatments, targeting both equal and unequal allocation. Initial treatment assignments (``j=1, 2, \ldots, m_0``) are made completely at random (``P_k(j)=\rho_k``, ``k=1, 2, \ldots, K``) until each group has at least one subject (i.e., ``N_k(m_0)>0``). Subsequent treatment assignments are made as follows:

```math
P_k(j) = \frac{\rho_k\left(\rho_k/\frac{N_k(j-1)}{j-1}\right)^\gamma}{\sum_{k=1}^K\rho_k\left(\rho_k/\frac{N_k(j-1)}{j-1}\right)^\gamma} ,\: j = m_0+1, \ldots, n.
```

where ``\gamma`` (``>0``) is a parameter of the procedure.

```@docs
DBCD
```

```@repl
using Incertus
dbcd = DBCD(2)     # DBCD, targeting 1:1 allocation (γ=2)

w = [1, 2, 3, 4]   # target allocation ratio
dbcd =  DBCD(w, 2) # DBCD, targeting 1:2:3:4 allocation (γ=2)
```

## Maximum Entropy Constraint Balance Randomization (MaxEnt)

A randomization procedure that can be used fo for two- and multi-arm trials with ``K`` treatments, targeting both equal and unequal allocation. 

Consider a point in the trial when ``j-1`` subjects have been randomized among the ``K``treatments, and let denote the corresponding treatment numbers as ``N_k(j-1)`` (``k=1, \ldots, K`` and ``\sum_{k_1}^KN_k(k-1) = j-1``). At the ``j^\text{th}`` allocation step, the randomization rule is defined as follows: 

a) For ``k=1, 2, \ldots, K``, compute ``B_k``, the hypothetical "lack of balnce", which results from assigning the ``j^\text{th}`` subject to treatment ``k``:

```math
B_k = \max\limits_{1\leq i \leq K}\left|\frac{N^{(k)}_i(j)}{k}-\rho_k\right|, \text{ where }N^{(k)}_i(j) = \left\{
\begin{array}{rl}
N_i(j-1) + 1, & i = k \\
N_i(j-1), & i \ne k
\end{array}
\right. .
```
b) The treatment randomization probabilities for the ``j^\text{th}`` subject ``\left(P_1(j), P_2(j), \ldots, P_K(j)\right)`` are determined as a solution to the constrained optimization problem:

```math
\begin{aligned}
&\text{minimize}\sum_{k=1}^KP_k(j)\log\left(\frac{P_k(j)}{\rho_k}\right) \\
&\text{subject to}\sum_{k=1}^KB_kP_k(j) \leq \eta B_{(1)} + (1-\eta)\sum_{k=1}^KB_k\rho_k \\
&\text{and}\sum_{k=1}^KP_k(j) = 1; \: 0\leq P_k(j) \leq 1, \: k = 1, 2, \ldots, K,
\end{aligned}
```

where ``B_{(1)}=\min\limits_{1\leq k \leq K}B_k``, and ``\eta`` ``\left(\in [0; 1]\right)`` is a parameter of the procedure that controls degree of randomness.

```@docs
MaxEnt
```

```@repl
using Incertus
maxent = MaxEnt(0.5)     # MaxEnt, targeting 1:1 allocation (η=2)

w = [1, 2, 3, 4]   # target allocation ratio
maxent = MaxEnt(w, 0.5)  # MaxEnt, targeting 1:2:3:4 allocation (η=2)
```


## Funcions implemented to calculate allocation probabilities

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

```@docs
allocation_prb(::MWUD, ::Vector{Int64})
```

```@docs
allocation_prb(::DBCD, ::Vector{Int64})
```

```@docs
allocation_prb(::MaxEnt, ::Vector{Int64})
```