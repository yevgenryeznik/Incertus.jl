function calc_trt_numbers(trt::Matrix{Int64})
    N = cumsum(trt, dims = 1)

    return(N)
end


function calc_imb(trt::Matrix{Int64})
    nsbj, ntrt = size(trt)

    @assert ntrt == 2 "The function is not implemented for multi-arm trials"
    
    # calculating treatment numbers 
    N = calc_trt_numbers(trt)

    # calculating imbalance
    imb = [N[i, 1] - N[i, 2] for i in 1:nsbj]

    return(imb)
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
    trt = sr.trt

    # getting simulation options
    nsbj, ntrt, nsim = size(trt)

    if ntrt == 2
        imb = [calc_imb(trt[:, :, s]) for s in 1:nsim]
        return [imb[s][nsbj] for s in 1:nsim]    
    else
        println("multi-arm trials are not supported yet")
        return nothing
    end
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
    trt = sr.trt
    ntrt = maximum(trt)
    if ntrt == 2
        trt = 2 .- trt
    end

    abs_imb = hcat([abs.(calc_imb(trt[:, :, s])) for s in axes(trt, 2)]...)
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
    trt = sr.trt
    ntrt = maximum(trt)
    if ntrt == 2
        trt = 2 .- trt
    end

    imb = hcat([calc_imb(trt[:, :, s]) for s in axes(trt, 2)]...)
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
    trt = sr.trt
    ntrt = maximum(trt)
    if ntrt == 2
        trt = 2 .- trt
    end

    abs_imb = hcat([abs.(calc_imb(trt[:, :, s])) for s in axes(trt, 2)]...)
    max_abs_imb = zeros(Int64, size(trt))
    for s in axes(trt, 2)
        for j in axes(trt, 1)
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