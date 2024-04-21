# Simulation example 1

Here, we simulate several randomization procedures, targeting 1:1 allocation.

The following procedures will be simulated:

- Comletele Randomized Design (`CRD`).
- Permuted Block Design with a block size equal to 2 (`PBD(1)`).
- Random Allocation Rule (`RAND`)
- Truncated Binomial Design (`TBD`)
- Big Stick Design with ``mti=3`` (BSD(3)).
- Efron's Biased Coin Design with ``p = \frac{2}{3}`` (`EBCD(2/3)`).
- Adjustable Biased Coin Design with ``a = 2`` (`ABCD(2)`).

We set sample size ``n=40`` and perform ``10,000`` simulations.

```@example
using Plots: savefig; nothing #hide
using Incertus

# sample size
nsbj = 40;

# number of simulations
nsim = 10000;

# randomization procedures to be simulated
rnd = [CRD(), PBD(1), RAND(nsbj), TBD(nsbj), BSD(3), EBCD(2//3), ABCD(2)];

# simulation run
sr = simulate(rnd, nsbj, nsim);

# calculating final imbalance, given simulations' output (`sr`) 
final_imb = calc_final_imb(sr);

# calculating expected absolute imbalance vs. allocation step, given simulations' output (`sr`) 
expected_abs_imb = calc_expected_abs_imb(sr);

# calculating variance of imbalance vs. allocation step, given simulations' output (`sr`) 
variance_of_imb = calc_variance_of_imb(sr);

# calculating expected maximum imbalances over first allocation steps, given simulations' output (`sr`) 
expected_max_abs_imb = calc_expected_max_abs_imb(sr);

# calculating cumulative average losses over first allocation steps, given simulations' output (`sr`) 
cummean_loss = calc_cummean_loss(sr);

# making a violin plot of final imbalances 
violin(final_imb)
savefig("violinplot.png"); nothing #hide

# making a plot of the expected absolute imbalances 
plot(expected_abs_imb, ylabel = "ecpexted absolute imbalance")
savefig("expected_absolute_imb_plot.png"); nothing #hide

# making a plot of the variances of imbalance
plot(variance_of_imb, ylabel = "variance of imbalance")
savefig("variance_of_imb_plot.png"); nothing #hide

# making a plot of expected maximum absolute imbalances over first allocation steps
plot(expected_max_abs_imb, ylabel = "expected maximum absolute imbalance")
savefig("expected_max_abs_imb_plot.png"); nothing #hide

# making a plot of cumulative average losses over the first allocation steps
plot(cummean_loss, ylabel = "cumulative average loss")
savefig("cummean_loss_plot.png"); nothing #hide
```

## Simulation results

### A _violin plot_ of the final imbalance after all treatment assignmnets complete

![](violinplot.png)

### Expected absolute imbalance vs. allocation step

![](expected_absolute_imb_plot.png)

### Variance of imbalance vs. allocation step

![](variance_of_imb_plot.png)

### Expected maximum absolute imbalance over first allocation steps

![](expected_max_abs_imb_plot.png)

### Cumulative average loss over first allocation steps

![](cummean_loss_plot.png)