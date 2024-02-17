# Simulation example

```@example
using Plots: savefig; nothing #hide
using Incertus

# sample size
nsbj = 40;

# number of simulations
nsim = 10000;

# randomization procedures to be simulated
rnd = [CRD(), PBD(1), BSD(3), EBCD(2//3)];

# simulation run
sr = simulate(rnd, nsbj, nsim);

# calculating final imbalance, given simulations' output (`sr`) 
final_imb = calc_final_imb(sr);

# making a violin plot of final imbalances 
violin(final_imb)
savefig("violinplot.png"); nothing #hide
```
!["Violin plot of the final imbalance after all treatment assignmnets complete."](violinplot.png)
