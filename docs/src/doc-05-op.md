# Operational characteristics

Several measures of imbalance and randomness have been implemented in the package:

## Measures of imbalance

Let us consider the following quantity as a measure of imbalance:

```math
D(j) = N_E(j)-N_C(j), \: \text{for two-arm trial and 1:1 randomization}, \: j = 1, 2, \ldots, n;
```

and

```math
D(j) = \sqrt{\sum\limits_{k=1}^K \left(N_k(j)-j\rho_k\right)^2}, \: \text{otherwise}, \: j = 1, 2, \ldots, n,
```

where the latter is the Euclidean distance between the _observed_ and the _targeted_ allocation.

The following measures of imbalance are available in the package.

### Final imbalance

``D(n)`` -- final _imbalance_/_distance (between the observed and targeted sample sizes)_ after all treatment assignments made.

```@docs
calc_final_imb(sr::SimulatedRandomization)
```

```@docs
calc_final_imb(sr::Vector{SimulatedRandomization})
```

### Expected absolute imbalance

``\mathbf{E}\left[|D(j)|\right]`` -- expected _absolute imbalance_/_distance (between the observed and targeted sample sizes)_ at the ``j^\text{th}`` allocation step.

```@docs
calc_expected_abs_imb(sr::SimulatedRandomization)
```

```@docs
calc_expected_abs_imb(sr::Vector{SimulatedRandomization})
```

### Variance of imbalance

``\mathbf{E}\left[|D(j)|^2\right]=\mathbf{var}\left[N_E(j)-N_C(j)\right]`` or ``=\sum_{k=1}^K\mathbf{var}\left[N_k(j)\right]`` -- _variance of imbalance_/_total variance of the treatment sample sizes_ at the ``j^\text{th}`` allocation step.

```@docs
calc_variance_of_imb(sr::SimulatedRandomization)
```

```@docs
calc_variance_of_imb(sr::Vector{SimulatedRandomization})
```

### Expected maximum imbalance over first allocation steps

``\mathbf{E}\left[\max\limits_{1\leq m \leq j}|D(m)|\right]`` -- expected maximum _imbalance_/_distance (between the observed and the targeted treatment sample sizes)_ over the first ``j`` allocation steps.

```@docs
calc_expected_max_abs_imb(sr::SimulatedRandomization)
```

```@docs
calc_expected_max_abs_imb(sr::Vector{SimulatedRandomization})
```

### A cumulative average loss over the first allocation steps

``Imb(j) = \frac{1}{j}\sum\limits_{m=1}^j\frac{\mathbf{E}\left[|D(m)|^2\right]}{m}`` -- a cumulative average loss at the ``j^\text{th}`` allocation step.

```@docs
calc_cummean_loss(sr::SimulatedRandomization)
```

```@docs
calc_cummean_loss(sr::Vector{SimulatedRandomization})
```

## Measures of randomness

### Cumulative average of expected proportions of correct guesses over first allocation steps under _two_ different guessing strategies 

``EPCG_{conv}(j)`` -- _cumulative average of expected proportions of correct guesses_ over first (``j``) allocation steps under the _convergence_ guessing strategy. 

- For a two-arm trial with ``1:1`` target allocation, it is defined as 
    
```math
\begin{aligned}
    EPCG_{conv}(j) &= \frac{1}{j}\sum\limits_{m=1}^j\mathbf{E}\left[G_m = \delta_m\right], \: j = 1, 2, \ldots, n, \\
    \text{where }G_m &= \left\{
        \begin{array}{rl}
            1, & D(m-1) < 0; \\
            \sim Bernoulli(0.5), & D(m-1) = 0; \\
            0, & D(m-1) > 0.
        \end{array}  
    \right. 
\end{aligned}    
```

- For a multu-arm trial, it is defined as 
    
```math
\begin{aligned}
    EPCG_{conv}(j) &= \frac{1}{j}\sum\limits_{m=1}^j\mathbf{E}\left[\mathbf{G}_m = \boldsymbol{\delta}_m\right], \: j = 1, 2, \ldots, n, \\
    \text{where }G_m &\sim Multinomial\left(1, \mathbf{P}\right), \\
    \text{where }\mathbf{P} &= \left(\frac{\left\{\Delta_1 = \min\limits_{1\leq i \leq K}\Delta_i\right\}}{\sum_{k=1}^K\left\{\Delta_k = \min\limits_{1\leq i \leq K}\Delta_i\right\}}, \ldots, \frac{\left\{\Delta_K = \min\limits_{1\leq i \leq K}\Delta_i\right\}}{\sum_{k=1}^K\left\{\Delta_k = \min\limits_{1\leq i \leq K}\Delta_i\right\}}\right), \\
    \text{where }\Delta_i &= \frac{N_i(m)}{m}-\rho_i, \: i = 1, 2, \ldots, K. 
\end{aligned}    
```

``G_m`` (or ``\mathbf{G}_m``) is a random variable taking values based on the investigator's guess; if ``G_m = \delta_m`` (``\mathbf{G}_m = \boldsymbol{\delta_m}``), the guess is correct at the ``m^\text{th}`` allocation step.

``EPCG_{max}(j)`` -- _cumulative average of expected proportions of correct guesses_ over first (``j``) allocation steps under the _maximum probability_ guessing strategy.

- For a two-arm trial with ``1:1`` target allocation, it is defined as 

```math
\begin{aligned}
EPCG_{max}(j) &= \frac{1}{j}\sum\limits_{m=1}^j\mathbf{E}\left[\widetilde{G}_m = \delta_m\right], \\
\text{where }\widetilde{G}_m &= \left\{
    \begin{array}{rl}
        1, & \phi_m > 0.5; \\
        \sim Bernoulli(0.5), & \phi_m = 0.5; \\
        0, & \phi_m < 0.5.
    \end{array}  
    \right. 
\end{aligned}
```

- For a multu-arm trial, it is defined as 
    
```math
\begin{aligned}
    EPCG_{max}(j) &= \frac{1}{j}\sum\limits_{m=1}^j\mathbf{E}\left[\widetilde{\mathbf{G}}_m = \boldsymbol{\delta}_m\right], \: j = 1, 2, \ldots, n, \\
    \text{where }G_m &\sim Multinomial\left(1, \mathbf{P}\right), \\
    \text{where }\mathbf{P} &= \left(\frac{\left\{P_1(m) = \max\limits_{1\leq i \leq K}P_i(m)\right\}}{\sum_{k=1}^K\left\{P_k(m) = \max\limits_{1\leq i \leq K}P_k(m)\right\}}, \ldots, \frac{\left\{P_K(m) = \max\limits_{1\leq i \leq K}P_i(m)\right\}}{\sum_{k=1}^K\left\{P_k(m) = \max\limits_{1\leq i \leq K}P_i(m)\right\}}\right). 
\end{aligned}    
```

``\widetilde{G}_m`` (or ``\widetilde{\mathbf{G}}_m``) is a random variable taking values based on the investigator's guess; if ``\widetilde{G}_m = \delta_m`` (``\widetilde{\mathbf{G}}_m = \boldsymbol{\delta_m}``), the guess is correct at the ``m^\text{th}`` allocation step.

```@docs
calc_cummean_epcg(sr::SimulatedRandomization, gs::String)
```

```@docs
calc_cummean_epcg(sr::Vector{SimulatedRandomization}, gs::String)
```

### Cumulative average of expected proportions of deterministic assignments over first allocation steps

``PD(j)`` -- _cumulative average of expected proportions of deterministic assignments_ over first (``j``) allocation steps.

- For a two-arm trial with ``1:1`` target allocation, it is defined as 

```math
PD(j) = \frac{1}{j}\sum\limits_{m=1}^j\Pr(\phi_m\in\{0, 1\})
```
- For a multu-arm trial, it is defined as 

```math
\begin{aligned}
PD(j) &= \frac{1}{j}\sum\limits_{m=1}^j\Pr(\mathbf{P}(m)\in\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_K\}), \\
\text{where }\mathbf{e}_k &= (0, \ldots, \underbrace{1}_{k^\text{th}\:\text{element}}, \ldots, 0), \: k = 1, 2, \ldots, K.
\end{aligned}
```

```@docs
calc_cummean_pda(sr::SimulatedRandomization)
```

```@docs
calc_cummean_pda(sr::Vector{SimulatedRandomization})
```

### Forcing index

``FI(j)`` -- _forcing index_, which takes values on a scale 0–1. 

- For a two-arm trial with ``1:1`` target allocation, it is defined as 

```math
FI(j) = \frac{4}{j}\sum\limits_{m=1}^j\mathbf{E}\left[|\phi_m-0.5|\right].
```

Note that ``FI(j) = 0, \forall j`` for CRD and ``FI(j) = 1`` for PBD with a block size ``bs=2``, assuming $j$ is even (most balanced design). 

- For a multu-arm trial, it is defined as 

```math
FI(j) = \frac{1}{j}\sum\limits_{m=1}^j\mathbf{E}\left[\sqrt{\sum\limits_{k=1}^K\left(P_k(m)-\rho_k\right)^2}\:\right].
```

Note that ``FI(j) = 0, \forall j`` for CRD. 

```@docs
calc_fi(sr::SimulatedRandomization)
```

```@docs
calc_fi(sr::Vector{SimulatedRandomization})
```


## Balance-randomness trade-off

``G(j) = \sqrt{\left\{Imb(j)\right\}^2 + \left\{FI(j)\right\}^2}`` represents a _balance-randomness trade-off_ at the ``j^\text{th}`` allocation step. Lower values of ``G(j)`` indicate better balance–randomness trade-off.

The following functions are available to deal with the characteristic.

```@docs
calc_brt(sr::SimulatedRandomization)
```

```@docs
calc_brt(sr::Vector{SimulatedRandomization})
```

## Allocation ratio preserving property

A procedure has an ARP property if

```math
\mathbf{E}\left[P_k(j)\right] = \rho_k, \: k = 1, 2, \ldots, K; \: j = 1, 2, \ldots, n,
```

where ``P_k(j)`` is the _conditional randomization probability_ for the ``k^\text{th}`` 
treatment group at the ``j^\text{th}`` allocation step, and ``\rho_k`` is a target allocation proportion for the ``k^\text{th}`` tretament.

The following functionality is avalable to evaluate ARP property:

```@docs
ARP
```

```@docs
eval_arp(sr::SimulatedRandomization)
```

```@docs
eval_arp(sr::Vector{SimulatedRandomization})
```