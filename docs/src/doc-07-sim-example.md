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

# calculating expected absolute imbalance vs. allocation step, given simulations' output (`sr`) 
expected_abs_imb = calc_expected_abs_imb(sr);

# making a violin plot of final imbalances 
violin(final_imb)
savefig("violinplot.png"); nothing #hide

# making a plot of the exåected absolute imbalances 
plot(expected_abs_imb)
savefig("expected_absolute_imb_plot.png"); nothing #hide
```

![Violin plot of the final imbalance after all treatment assignmnets complete.](violinplot.png)

![expected absolute imbalance vs. allocation step](expected_absolute_imb_plot.png)
