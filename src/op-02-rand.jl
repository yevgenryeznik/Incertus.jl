# calculating correct guess probability under the convergence guessing strategy
function calc_guess_cgs(trt::Matrix{Int64}, target::Vector{<:Number})
    nsbj, ntrt = size(trt)
    
    # calculating treatment numbers
    N = calc_trt_numbers(trt)

    if (ntrt == 2) && (target[1] == target[2])
        # calculating imbalance before every allocation step
        imb = [0; N[1:end-1, 1] - N[1:end-1, 2]]

        # making guess, given current imbalance
        guess = [d == 0 ? rand(Binomial(1, 0.5)) : (d < 0 ? 1 : 0) for d in imb]

        # returning correct guesses
        return Int.(guess .== trt[:, 1])
    else
        # target allocation proportions
        ρ = target ./ sum(target)

        # differences between current and target allocation proportions
        Δ = [zeros(Int64, 1, ntrt); N[1:end-1, :] ./ (1:nsbj-1) .- ρ']

        # making guess, given current imbalance
        guess = zeros(Int64, size(trt))
        for j in axes(guess, 1)
            # guessing probabilities, given treatment numbers
            pnum = Int.(Δ[j, :] .== minimum(Δ[j, :]))
            p = pnum ./ sum(pnum)

            # making a guess
            guess[j, :] = rand(Multinomial(1, p))
        end

        # returning correct guesses
        return [Int.(guess[j, :] == trt[j, :]) for j in 1:nsbj]
    end
end


# calculating correct guess probability under the maximum probability guessing strategy
function calc_guess_mpgs(trt::Matrix{Int64}, prb::Matrix{Float64}, target::Vector{<:Number})
    nsbj, ntrt = size(trt)
    
    if (ntrt == 2) && (target[1] == target[2])
        # making guess, given allocation probability
        guess = [p == 0.5 ? rand(Binomial(1, 0.5)) : (p > 0.5 ? 1 : 0) for p in prb[:, 1]]

        # returning correct guesses
        return Int.(guess .== trt[:, 1])
    else
        # making guess, given current allocation probabilities
        guess = zeros(Int64, size(trt))
        for j in axes(guess, 1)
            # guessing probabilities, given treatment numbers
            pnum = Int.(prb[j, :] .== maximum(prb[j, :]))
            p = pnum ./ sum(pnum)

            # making a guess
            guess[j, :] = rand(Multinomial(1, p))
        end

        # returning correct guesses
        return [Int.(guess[j, :] == trt[j, :]) for j in 1:nsbj]
    end
end


# calculating deterministic assignments
function calc_da(prb::Matrix{Float64}, target::Vector{<:Number})
    nsbj, ntrt = size(trt)
    
    if (ntrt == 2) && (target[1] == target[2])
        # returning deterministic assignment 
        return [Int(prb[j, 1] ∈ [0, 1]) for j in 1:nsbj]
    else
        # returning deterministic assignment 
        return [Int(any(prb[j, :] .≈ 1)) for j in 1:nsbj]
    end
end


"""Function calculates cumulative averages of expected proportions of correct guesses vs. allocation step.

# Call
`calc_cummean_epcg(sr, gs)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.
- `gs::String`: guessing strategy; accepts two values: `"C"` (corresponds to the _convergence_ guessing strategy) of `"MP"` (corresponds to the _maximum probability_ guessing strategy).

# Result
- A `Vector` of _cumulative averages of the expected proportions of correct guesses_ summarized via simulations.
"""
function calc_cummean_epcg(sr::SimulatedRandomization, gs::String)
    @assert gs in ["C", "MP"] "`gs` input parameter must have one of the following values: \"C\" or \"MP\"."
    
    # target allocation ratio
    target = sr.target

    # extracting treatment assignments
    trt = sr.trt

    # extracting probabilities of treatment assignments
    prb = sr.prb
    
    nsbj, _, nsim = size(trt)
    
    pcg = gs == "C" ? 
        hcat([calc_guess_cgs(trt[:, :, s], target) for s in 1:nsim]...) :
        hcat([calc_guess_mpgs(trt[:, :, s], prb[:, :, s], target) for s in 1:nsim]...)
    
    expected_pcg = vec(mean(pcg, dims = 2))
    sbj = collect(1:nsbj)
    cummean_epcg = cumsum(expected_pcg) ./ sbj

    return cummean_epcg
end


"""Function calculates cumulative averages of expected proportions of correct guesses vs. allocation step.

# Call
`calc_cummean_epcg(sr, gs)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.
- `gs::String`: guessing strategy; accepts two values: `"C"` (corresponds to the _convergence_ guessing strategy) of `"MP"` (corresponds to the _maximum probability_ guessing strategy).

# Result
- A vector of _cumulative averages of the expected proportions of correct guesses_ summarized via simulations.
"""
function calc_cummean_epcg(sr::Vector{SimulatedRandomization}, gs::String)
    cummean_epcg = hcat([calc_cummean_epcg(item, gs) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(cummean_epcg, rnd_labels)
end



"""Function calculates cumulative averages of the proportions of deterministic assignments vs. allocation step.

# Call
`calc_cummean_pda(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.

# Result
- A vector of the _cumulative averages of the proportions of deterministic assignmnets_ values summarized via simulations.
"""
function calc_cummean_pda(sr::SimulatedRandomization)
    # target allocation ratio
    target = sr.target

    # extracting treatment assignments
    trt = sr.trt

    # extracting probabilities of treatment assignments
    prb = sr.prb
    
    nsbj, _, nsim = size(trt)

    da = hcat([calc_da(prb[:, :, s], target) for s in 1:nsim]...)
    
    pda = vec(mean(da, dims = 2))
    sbj = collect(1:nsbj)
    cummean_pda = cumsum(pda) ./ sbj

    return cummean_pda
end


"""Function calculates cumulative averages of the proportions of deterministic assignments vs. allocation step.

# Call
`calc_cummean_pda(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of the _cumulative averages of the proportions of deterministic assignmnets_ summarized via simulations.
"""
function calc_cummean_pda(sr::Vector{SimulatedRandomization})
    cummean_pda = hcat([calc_cummean_pda(item) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(cummean_pda, rnd_labels)
end


"""Function calculates forcing index vs. allocation step.

# Call
`calc_fi(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.

# Result
- A `Vector` of the _forcing index_ values summarized via simulations.
"""
function calc_fi(sr::SimulatedRandomization)
    # target allocation ratio
    target = sr.target

    # extracting treatment assignments
    trt = sr.trt

    # extracting probabilities of treatment assignments
    prb = sr.prb
    
    nsbj, ntrt, nsim = size(trt)

    if (ntrt == 2) && (target[1] == target[2])
        fi1 = hcat([abs.(prb[:, 1, s] .- 0.5) for s in 1:nsim]...)
        efi1 = vec(mean(fi1, dims = 2))
        sbj = collect(1:nsbj)
        
        return 4 .* (cumsum(efi1) ./ sbj)
    else
        # target allocation proportion
        ρ = target ./ sum(target)

        fi1 = hcat([sqrt((prb[:, :, s] - ρ)'*(prb[:, :, s] - ρ)) for s in 1:nsim]...)
        efi1 = vec(mean(fi1, dims = 2))
        sbj = collect(1:nsbj)
        
        return cumsum(efi1) ./ sbj
    end
end


"""Function calculates forcing index vs. allocation step.

# Call
`calc_fi(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of the _forcing index_ values summarized via simulations.
"""
function calc_fi(sr::Vector{SimulatedRandomization})
    fi = hcat([calc_fi(item) for item in sr]...)
    rnd_labels = [item.label for item in sr]

    return DataFrame(fi, rnd_labels)
end