function calc_trt_numbers(trt::Matrix{Int64})
    N = cumsum(trt, dims = 1)

    return(N)
end


function calc_imb(trt::Matrix{Int64}, target::Vector{<:Number})
    nsbj, ntrt = size(trt)

    # calculating treatment numbers 
    N = calc_trt_numbers(trt)

    if (ntrt == 2) & (target[1] == target[2])
        # calculating imbalance
        return [N[j, 1] - N[j, 2] for j in 1:nsbj]
    else
        # target allocation proportions
        ρ = target ./ sum(target)

        # calculating (Euclidean) distance
        return [sqrt((N[j, :] - j .* ρ)'*(N[j, :] - j .* ρ)) for j in 1:nsbj]
    end
end


"""Function calculates final imbalance

# Call
`calc_final_imb(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.

# Result
- A `Vector` of _final imbalance_ values obtained via simulations.
"""
function calc_final_imb(sr::SimulatedRandomization)
    # target allocation ratio
    target = sr.target

    # extracting treatment assignments
    trt = sr.trt

    # getting simulation options
    _, _, nsim = size(trt)

    imb = [calc_imb(trt[:, :, s], target) for s in 1:nsim]
    return [imb[s][end] for s in 1:nsim]    
end


"""Function calculates final imbalance

# Call
`calc_final_imb(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of simulated _final imbalances_' values obtained via simulations.
"""
function calc_final_imb(sr::Vector{SimulatedRandomization})
    final_imb = hcat([calc_final_imb(item) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(final_imb, rnd_labels)
end


"""Function calculates expected absolute imbalance vs. allocation step.

# Call
`calc_expected_abs_imb(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.

# Result
- A `Vector` of _expected absolute imbalance_ values summarized via simulations.
"""
function calc_expected_abs_imb(sr::SimulatedRandomization)
    # target allocation ratio
    target = sr.target
    
    # extracting treatment assignments
    trt = sr.trt

    # getting simulation options
    _, _, nsim = size(trt)

    abs_imb = hcat([abs.(calc_imb(trt[:, :, s], target)) for s in 1:nsim]...)
    expected_abs_imb = mean(abs_imb, dims = 2)

    return vec(expected_abs_imb)
end


"""Function calculates expected absolute imbalance vs. allocation step.

# Call
`calc_expected_abs_imb(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of _expected absolute imbalances_' values summarized via simulations.
"""
function calc_expected_abs_imb(sr::Vector{SimulatedRandomization})
    expected_abs_imb = hcat([calc_expected_abs_imb(item) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(expected_abs_imb, rnd_labels)
end


"""Function calculates variance of imbalance vs. allocation step.

# Call
`calc_variance_of_imb(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.

# Result
- A `Vector` of _variance of imbalance_ values summarized via simulations.
"""
function calc_variance_of_imb(sr::SimulatedRandomization)
    # target allocation ratio
    target = sr.target
    
    # extracting treatment assignments
    trt = sr.trt

    # getting simulation options
    _, _, nsim = size(trt)

    imb = hcat([calc_imb(trt[:, :, s], target) for s in 1:nsim]...)
    imb2 = imb.^2
    variance_of_imb = mean(imb2, dims = 2)

    return vec(variance_of_imb)
end


"""Function calculates variance of imbalance vs. allocation step.

# Call
`calc_variance_of_imb(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of _variances' of imbalance_ values summarized via simulations.
"""
function calc_variance_of_imb(sr::Vector{SimulatedRandomization})
    variance_of_imb = hcat([calc_variance_of_imb(item) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(variance_of_imb, rnd_labels)
end


"""Function calculates expected maximum absolute imbalance over first allocations vs. allocation step.

# Call
`calc_expected_max_abs_imb(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.

# Result
- A `Vector` of _expected maximum absolute imbalance over firat allocations_ values summarized via simulations.
"""
function calc_expected_max_abs_imb(sr::SimulatedRandomization)
    # target allocation ratio
    target = sr.target
    
    # extracting treatment assignments
    trt = sr.trt

    # getting simulation options
    nsbj, _, nsim = size(trt)

    abs_imb = hcat([abs.(calc_imb(trt[:, :, s], target)) for s in 1:nsim]...)
    max_abs_imb = zeros(Number, nsbj, nsim)
    for s in 1:nsim
        for j in 1:nsbj
            max_abs_imb[j, s] = maximum(abs_imb[1:j, s])
        end
    end
    expected_max_abs_imb = mean(max_abs_imb, dims = 2)

    return vec(expected_max_abs_imb)
end


"""Function calculates expected maximum absolute imbalance over first allocations vs. allocation step.

# Call
`calc_expected_max_abs_imb(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of _expected maximum absolute imbalances'_ values (over first allocations) summarized via simulations.
"""
function calc_expected_max_abs_imb(sr::Vector{SimulatedRandomization})
    expected_max_abs_imb = hcat([calc_expected_max_abs_imb(item) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(expected_max_abs_imb, rnd_labels)
end


"""Function calculates cumulative average loss vs. allocation step.

# Call
`calc_cummean_loss(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.

# Result
- A `Vector` of _cumulative average loss'_ values summarized via simulations.
"""
function calc_cummean_loss(sr::SimulatedRandomization)
    
    variance_of_imb = calc_variance_of_imb(sr)
    sbj = eachindex(variance_of_imb)

    cummean_loss = cumsum(variance_of_imb ./ sbj) ./sbj
    return cummean_loss
end


"""Function calculates cumulative average loss vs. allocation step.

# Call
`calc_cummean_loss(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of _cumulative average losses'_ values summarized via simulations.
"""
function calc_cummean_loss(sr::Vector{SimulatedRandomization})
    cummean_loss = hcat([calc_cummean_loss(item) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(cummean_loss, rnd_labels)
end