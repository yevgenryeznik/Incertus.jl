function calc_imb(δ::Vector{Int64})
    imb = 2 .* cumsum(δ) - eachindex(δ)

    return(imb)
end


"""Function used to calculate expected absolute imbalance vs. allocation step.

# Call
`calc_expected_abs_imb(sr)`

# Arguments
- sr::SimulatedRandomization: an instance of `SimulatedRandomization`, an object, 
representing simulation output.

# Result
- A vector of _expected absolute imbalance_ values summarized via simulations.
"""
function calc_expected_abs_imb(sr::SimulatedRandomization)
    trt = sr.trt
    ntrt = maximum(2)
    if ntrt == 2
        trt = 2 .- trt
    end

    abs_imb = hcat([abs.(calc_imb(trt[:, s])) for s in axes(trt, 2)]...)
    expected_abs_imb = mean(abs_imb, dims = 2)

    return vec(expected_abs_imb)
end

function calc_variance_of_imb(sr::SimulatedRandomization)
    trt = sr.trt
    ntrt = maximum(2)
    if ntrt == 2
        trt = 2 .- trt
    end

    imb = hcat([calc_imb(trt[:, s]) for s in axes(trt, 2)]...)
    imb2 = imb.^2
    variance_of_imb = mean(imb2, dims = 2)

    return vec(variance_of_imb)
end

function calc_expected_max_abs_imb(sr::SimulatedRandomization)
    trt = sr.trt
    ntrt = maximum(2)
    if ntrt == 2
        trt = 2 .- trt
    end

    abs_imb = hcat([abs.(calc_imb(trt[:, s])) for s in axes(trt, 2)]...)
    max_abs_imb = zeros(Int64, size(trt))
    for s in axes(trt, 2)
        for j in axes(trt, 1)
            max_abs_imb[j, s] = maximum(abs_imb[1:j, s])
        end
    end
    expected_max_abs_imb = mean(max_abs_imb, dims = 2)

    return vec(expected_max_abs_imb)
end

function calc_cummean_loss(sr::SimulatedRandomization)
    
    variance_of_imb = calc_variance_of_imb(sr)
    sbj = eachindex(variance_of_imb)

    cummean_loss = cumsum(variance_of_imb ./ sbj) ./sbj
    return cummean_loss
end