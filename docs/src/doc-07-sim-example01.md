# Simulation example 1

## Simulation setup

Here, we simulate several randomization procedures, targeting 1:1 allocation.

The following procedures will be simulated:

- Comletele Randomized Design (`CRD`).
- Permuted Block Design with a block size equal to 2 (`PBD(1)`).
- Random Allocation Rule (`RAND`).
- Truncated Binomial Design (`TBD`).
- Big Stick Design with ``mti=3`` (`BSD(3)`).
- Efron's Biased Coin Design with ``p = \frac{2}{3}`` (`EBCD(2/3)`).
- Adjustable Biased Coin Design with ``a = 2`` (`ABCD(2)`).

We set sample size ``n=40`` and perform ``10,000`` simulations.

## Simulation run

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

# ========== Calculating operational characteristics ==========
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

# cumulative average of expected proportions of correct guesses over first allocation steps under the convergence guessing strategy, 
# given simulations' output (`sr`) 
cummean_epcg_c = calc_cummean_epcg(sr, "C");

# cumulative average of expected proportions of correct guesses over first allocation steps under the maximum probability guessing strategy, 
# given simulations' output (`sr`) 
cummean_epcg_mp = calc_cummean_epcg(sr, "MP");

# cumulative average of expected proportions of deterministic assignments over first alocation steps, given simulations' output (`sr`) 
cummean_pda = calc_cummean_pda(sr);

# calculating forcing index vs. allocation step, given simulations' output (`sr`) 
fi = calc_fi(sr);

# Evaluating balance-randomness tradeoff vs. allocation step
brt = calc_brt(sr);


# ========== Visualizing results ==========
if !isdir("./figures01") # hide
    mkdir("./figures01") # hide
end                      # hide
# making a violin plot of final imbalances 
violin(final_imb)
savefig("./figures01/01-final-imb-violin-plot.png"); nothing #hide

# making a plot of the expected absolute imbalances 
plot(expected_abs_imb, ylabel = "ecpexted absolute imbalance")
savefig("./figures01/02-expected-absolute-imb-plot.png"); nothing #hide

# making a plot of the variances of imbalance
plot(variance_of_imb, ylabel = "variance of imbalance")
savefig("./figures01/03-variance-of-imb-plot.png"); nothing #hide

# making a plot of expected maximum absolute imbalances over first allocation steps
plot(expected_max_abs_imb, ylabel = "expected maximum absolute imbalance")
savefig("./figures01/04-expected-max-abs-imb-plot.png"); nothing #hide

# making a plot of cumulative average losses over the first allocation steps
plot(cummean_loss, ylabel = "cumulative average loss")
savefig("./figures01/05-cummean-loss-plot.png"); nothing #hide

# making a plot of cumulative averages of expected proportions of correct guesses over first allocation steps 
# (under the convergence guessing strategy) 
plot(cummean_epcg_c, ylabel = "cumulative average of EPCGs (CGS)")
savefig("./figures01/06-cummean-epcg-c-plot.png"); nothing #hide

# making a plot of cumulative averages of expected proportions of correct guesses over first allocation steps 
# (under the maximum probability guessing strategy) 
plot(cummean_epcg_mp, ylabel = "cumulative average of EPCGs (MPGS)")
savefig("./figures01/07-cummean-epcg-mp-plot.png"); nothing #hide

# making a plot of cumulative averages of expected proportions of deterministic assignments over first alocation steps
plot(cummean_pda, ylabel = "cumulative average of PDAs")
savefig("./figures01/08-cummean-pda-plot.png"); nothing #hide

# making a plot of forcing indeces vs. allocation step 
plot(fi, ylabel = "forcing index")
savefig("./figures01/09-fi-plot.png"); nothing #hide

# making a plot of forcing indeces vs. allocation step 
heatmap(brt)
savefig("./figures01/10-brt-heatmap-plot.png"); nothing #hide

```

## Simulation results

### A _violin plot_ of the final imbalance after all treatment assignmnets complete

![](./figures01/01-final-imb-violin-plot.png)

### Expected absolute imbalance vs. allocation step

![](./figures01/02-expected-absolute-imb-plot.png)

### Variance of imbalance vs. allocation step

![](./figures01/03-variance-of-imb-plot.png)

### Expected maximum absolute imbalance over first allocation steps

![](./figures01/04-expected-max-abs-imb-plot.png)

### Cumulative average loss over first allocation steps

![](./figures01/05-cummean-loss-plot.png)

### Cumulative averages of expected proportions of correct guesses over first allocation steps (under the convergence guessing strategy) 

![](./figures01/06-cummean-epcg-c-plot.png)

### Cumulative averages of expected proportions of correct guesses over first allocation steps (under the maximum probability guessing strategy) 

![](./figures01/07-cummean-epcg-mp-plot.png)

### Cumulative averages of expected proportions of deterministic assignments over first alocation steps

![](./figures01/08-cummean-pda-plot.png)

### Forcing indeces vs. allocation step 

![](./figures01/09-fi-plot.png)

### Balance-randomness tradeoff vs. allocation step 

![](./figures01/10-brt-heatmap-plot.png)
