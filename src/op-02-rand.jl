# calculating correct guess probability under the convergence guessing strategy
function calc_guess_cgs(δ::Vector{Int64})
    # calculating imbalance before every allocation step
    imb = [0; 2 .* cumsum(δ[1:end-1]) - eachindex(δ[1:end-1])]

    # making guess, given current imbalance
    guess = [d == 0 ? rand(Binomial(1, 0.5)) : (d < 0 ? 1 : 0) for d in imb]

    # correct guesses
    correct_guess = Int.(guess .== δ)

    return(correct_guess)
end


# calculating correct guess probability under the maximum probability guessing strategy
function calc_guess_mpgs(δ::Vector{Int64}, ϕ::Vector{Float64})
    # making guess, given allocation probability
    guess = [p == 0.5 ? rand(Binomial(1, 0.5)) : (p > 0.5 ? 1 : 0) for p in ϕ]

    # correct guesses
    correct_guess = Int.(guess .== δ)

    return(correct_guess)
end


# calculating deterministic assignments
function calc_da(ϕ::Vector{Float64})
    deterministic_assignment = [Int(item ∈ [0, 1]) for item in ϕ]

    return(deterministic_assignment)
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
    
    trt = sr.trt
    prb = sr.prb
    
    nsbj, ntrt, nsim = size(trt)
    
    if ntrt == 2
        pcg = gs == "C" ? 
            hcat([calc_guess_cgs(trt[:, 1, s]) for s in 1:nsim]...) :
            hcat([calc_guess_mpgs(trt[:, 1, s], prb[:, 1, s]) for s in 1:nsim]...)
    
        expected_pcg = vec(mean(pcg, dims = 2))
        sbj = collect(1:nsbj)
        cummean_epcg = cumsum(expected_pcg) ./ sbj

        return cummean_epcg
    else
        println("multi-arm trials are not supported yet")
        return nothing
    end
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
    trt = sr.trt
    prb = sr.prb
    
    nsbj, ntrt, nsim = size(trt)

    if ntrt == 2
        da = hcat([calc_da(prb[:, 1, s]) for s in 1:nsim]...)
    
        pda = vec(mean(da, dims = 2))
        sbj = collect(1:nsbj)
        cummean_pda = cumsum(pda) ./ sbj

        return cummean_pda
    else
        println("multi-arm trials are not supported yet")
        return nothing
    end
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
    trt = sr.trt
    prb = sr.prb
    
    nsbj, ntrt, nsim = size(trt)

    if ntrt == 2
        fi1 = hcat([abs.(prb[:, 1, s] .- 0.5) for s in 1:nsim]...)
    
        efi1 = vec(mean(fi1, dims = 2))
        sbj = collect(1:nsbj)
        fi = 4 .* (cumsum(efi1) ./ sbj)

        return fi
    else
        println("multi-arm trials are not supported yet")
        return nothing    
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