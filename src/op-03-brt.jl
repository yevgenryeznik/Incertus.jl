"""Function calculates balance-randomness trade-off vs. allocation step.

# Call
`calc_brt(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, 
representing simulation output.

# Result
- A vector of _balance-randomness trade-off measurement_ values summarized via simulations.
"""
function calc_brt(sr::SimulatedRandomization)
    # measure of imbalance vs. allocation step
    cummean_loss = calc_cummean_loss(sr)

    # measure of randomness vs. allocation step
    fi = calc_fi(sr)

    # balance-randomness trade-off vs. allocation step
    G = [sqrt(item1^2 + item2^2) for (item1, item2) in zip(cummean_loss, fi)]

    return G
end