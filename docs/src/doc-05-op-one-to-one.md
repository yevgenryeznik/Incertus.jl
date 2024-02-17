# Operational characteristics (_two-arm trial with equal 1:1 allocation_)

Several measures of imbalance and randomness have been implemented in the package:

_**Measures of imbalance**_:

- ``D(n)`` -- final imbalance after all treatment assignments made.
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

```@docs
calc_variance_of_imb(sr::Vector{SimulatedRandomization})
```

## Expected maximum imbalance over first allocation steps

```@docs
calc_expected_max_abs_imb(sr::SimulatedRandomization)
```

```@docs
calc_expected_max_abs_imb(sr::Vector{SimulatedRandomization})
```

## A cumulative average loss at the first allocation steps

```@docs
calc_cummean_loss(sr::SimulatedRandomization)
```

```@docs
calc_cummean_loss(sr::Vector{SimulatedRandomization})
```

## Cumulative average of expected proportions of correct guesses over first allocation steps under _two_ different guessing strategies 

```@docs
calc_cummean_epcg(sr::SimulatedRandomization, gs::String)
```

```@docs
calc_cummean_epcg(sr::Vector{SimulatedRandomization}, gs::String)
```

## Cumulative average of expected proportions of deterministic assignments over first allocation steps

```@docs
calc_cummean_pda(sr::SimulatedRandomization)
```

```@docs
calc_cummean_pda(sr::Vector{SimulatedRandomization})
```

## Forcing index

```@docs
calc_fi(sr::SimulatedRandomization)
```

```@docs
calc_fi(sr::Vector{SimulatedRandomization})
```

## Balance-randomness trade-off

```@docs
calc_brt(sr::SimulatedRandomization)
```

```@docs
calc_brt(sr::Vector{SimulatedRandomization})
```