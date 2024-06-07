# Simulation example 2

## Simulation setup

Here, we simulate several randomization procedures, targeting 4:3:2:1 allocation.

The following procedures will be simulated:

DLUD(``a=2``), and MWUD(``\alpha=2``)
- Comletele Randomized Design (`CRD`).
- Permuted Block Design with a block size equal to 2 (`PBD(1)`).
- Block Urn Design with parameter ``\lambda=2`` (`BUD(2)`).
- Random Allocation Rule (`RAND`).
- Truncated Multinomial Design (`TMD`).
- Drop-the-Loser Urn Design with ``a=2`` (`DLUD(a=2)`.
- Mass Weighted Urn Design with ``\alpha=2`` (`MWUD(2)`).

We set sample size ``n=40`` and perform ``10,000`` simulations.

## Simulation run

```@example
using Plots: savefig; nothing #hide
using Incertus

# sample size
nsbj = 40;

# number of simulations
nsim = 10000;

# target allocation
w = [4, 3, 2, 1];

# randomization procedures to be simulated
rnd = [CRD(w), PBD(w, 1), BUD(w, 2), RAND(w, nsbj), TMD(w, nsbj), DLUD(w, 2), MWUD(w, 2)];

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

# evaluating allocation ratio preserving (ARP) property 
arp = eval_arp(sr);

# ========== Visualizing results ==========
if !isdir("./figures02") # hide
    mkdir("./figures02") # hide
end                      # hide
# making a violin plot of final imbalances 
violin(final_imb)
savefig("./figures02/01-final-imb-violin-plot.png"); nothing #hide

# making a plot of the expected absolute imbalances 
plot(expected_abs_imb, ylabel = "ecpexted absolute imbalance")
savefig("./figures02/02-expected-absolute-imb-plot.png"); nothing #hide

# making a plot of the variances of imbalance
plot(variance_of_imb, ylabel = "variance of imbalance")
savefig("./figures02/03-variance-of-imb-plot.png"); nothing #hide

# making a plot of expected maximum absolute imbalances over first allocation steps
plot(expected_max_abs_imb, ylabel = "expected maximum absolute imbalance")
savefig("./figures02/04-expected-max-abs-imb-plot.png"); nothing #hide

# making a plot of cumulative average losses over the first allocation steps
plot(cummean_loss, ylabel = "cumulative average loss")
savefig("./figures02/05-cummean-loss-plot.png"); nothing #hide

# making a plot of cumulative averages of expected proportions of correct guesses over first allocation steps 
# (under the convergence guessing strategy) 
plot(cummean_epcg_c, ylabel = "cumulative average of EPCGs (CGS)")
savefig("./figures02/06-cummean-epcg-c-plot.png"); nothing #hide

# making a plot of cumulative averages of expected proportions of correct guesses over first allocation steps 
# (under the maximum probability guessing strategy) 
plot(cummean_epcg_mp, ylabel = "cumulative average of EPCGs (MPGS)")
savefig("./figures02/07-cummean-epcg-mp-plot.png"); nothing #hide

# making a plot of cumulative averages of expected proportions of deterministic assignments over first alocation steps
plot(cummean_pda, ylabel = "cumulative average of PDAs")
savefig("./figures02/08-cummean-pda-plot.png"); nothing #hide

# making a plot of forcing indeces vs. allocation step 
plot(fi, ylabel = "forcing index")
savefig("./figures02/09-fi-plot.png"); nothing #hide

# making a plot of forcing indeces vs. allocation step 
heatmap(brt)
savefig("./figures02/10-brt-heatmap-plot.png"); nothing #hide

# making a plot to demonstrate ARP property
plot(arp)
savefig("./figures02/11-arp-plot.png"); nothing #hide

nothing
```

## Simulation results

### A _violin plot_ of the final imbalance after all treatment assignmnets complete

![](./figures02/01-final-imb-violin-plot.png)

### Expected absolute imbalance vs. allocation step

![](./figures02/02-expected-absolute-imb-plot.png)

### Variance of imbalance vs. allocation step

![](./figures02/03-variance-of-imb-plot.png)

### Expected maximum absolute imbalance over first allocation steps

![](./figures02/04-expected-max-abs-imb-plot.png)

### Cumulative average loss over first allocation steps

![](./figures02/05-cummean-loss-plot.png)

### Cumulative averages of expected proportions of correct guesses over first allocation steps (under the convergence guessing strategy) 

![](./figures02/06-cummean-epcg-c-plot.png)

### Cumulative averages of expected proportions of correct guesses over first allocation steps (under the maximum probability guessing strategy) 

![](./figures02/07-cummean-epcg-mp-plot.png)

### Cumulative averages of expected proportions of deterministic assignments over first alocation steps

![](./figures02/08-cummean-pda-plot.png)

### Forcing indeces vs. allocation step 

![](./figures02/09-fi-plot.png)

### Balance-randomness tradeoff vs. allocation step 

![](./figures02/10-brt-heatmap-plot.png)

### Allocation ratio preserving (ARP) property vs. allocation step 

![](./figures02/11-arp-plot.png)
