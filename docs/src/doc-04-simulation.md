# Simulation

To perform simulations, the following functionality has been implemented:

```@docs
SimulatedRandomization
```

```@docs
simulate(rnd::T, nsbj::Int64, nsim::Int64, seed::Int64 = 314159) where T <: Randomization
```

```@docs
simulate(rnd::Vector{<:Randomization}, nsbj::Int64, nsim::Int64, seed::Int64 = 314159)
```