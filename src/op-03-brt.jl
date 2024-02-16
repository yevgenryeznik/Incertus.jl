"""Function calculates balance-randomness trade-off vs. allocation step.

# Call
`calc_brt(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, 
representing simulation output.

# Result
- A `Vector` of _balance-randomness trade-off measurement_ values summarized via simulations.
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


"""Function calculates balance-randomness trade-off vs. allocation step.

# Call
`calc_brt(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of _balance-randomness trade-off measurements_' values summarized via simulations.
"""
function calc_brt(sr::Vector{SimulatedRandomization})
    brt = hcat([calc_brt(item) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(brt, rnd_labels)
end