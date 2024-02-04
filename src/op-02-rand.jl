# calculating correct guess probability under the convergence guessing strategy
function calc_guess_cgs(δ::Vector{Int64})
    imb = [0; 2 .* cumsum(δ[1:end-1]) - eachindex(δ[1:end-1])]
    G = [0.5*(1-sign(item)) for item in imb]
    guess = [rand(Binomial(1, p)) for p in G]

    return(Int.(guess .== δ))
end


# calculating correct guess probability under the maximum probability guessing strategy
function calc_guess_mpgs(ϕ::Vector{Float64})
    guess = [0.5*(1-sign(item-0.5)) for item in ϕ]

    return(guess)
end


# calculating deterministic assignments
function calc_da(ϕ::Vector{Float64})
    da = [Int(item ∈ [0, 1]) for item in ϕ]

    return(da)
end

"""Function calculates cumulative averages of expected proportions of correct guesses vs. allocation step.

# Call
`calc_cummean_epcg(sr, gs)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, 
representing simulation output.
- `gs::String`: guessing strategy; accepts two values: `"C"` (corresponds to the 
_convergence_ guessing strategy) of `"MP"` (corresponds to the _maximum probability_ 
guessing strategy).

# Result
- A vector of _cumulative averages of the expected proportions of correct guesses_ values summarized via simulations.
"""
function calc_cummean_epcg(sr::SimulatedRandomization, gs::String)
    @assert gs in ["C", "MP"] "`gs` input parameter must have of the values in [\"C\", \"MP\"]"
    trt = sr.trt
    prb = sr.prb
    ntrt = maximum(trt)
    if ntrt == 2
        trt = 2 .- trt
        prb = hcat([prb[:, 1, s] for s in axes(prb, 3)]...) 
    end

    pcg = gs == "C" ? 
        hcat([calc_guess_cgs(trt[:, s]) for s in axes(trt, 2)]...) :
        hcat([calc_guess_mpgs(prb[:, s]) for s in axes(trt, 2)]...)
    
    expected_pcg = vec(mean(pcg, dims = 2))
    sbj = eachindex(expected_pcg)
    cummean_epcg = cumsum(expected_pcg) ./ sbj

    return cummean_epcg
end


"""Function calculates cumulative averages of the proportions of deterministic assignments vs. allocation step.

# Call
`calc_cummean_pda(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, 
representing simulation output.

# Result
- A vector of the _cumulative averages of the proportions of deterministic assignmnets_ values summarized via simulations.
"""
function calc_cummean_pda(sr::SimulatedRandomization)
    trt = sr.trt
    prb = sr.prb
    ntrt = maximum(trt)
    if ntrt == 2
        trt = 2 .- trt
        prb = hcat([prb[:, 1, s] for s in axes(prb, 3)]...) 
    end

    da = hcat([calc_da(prb[:, s]) for s in axes(prb, 2)]...)
    
    pda = vec(mean(da, dims = 2))
    sbj = eachindex(pda)
    cummean_pda = cumsum(pda) ./ sbj

    return cummean_pda
end


"""Function calculates forcing index vs. allocation step.

# Call
`calc_fi(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, 
representing simulation output.

# Result
- A vector of the _forcing index_ values summarized via simulations.
"""
function calc_fi(sr::SimulatedRandomization)
    trt = sr.trt
    prb = sr.prb
    ntrt = maximum(trt)
    if ntrt == 2
        trt = 2 .- trt
        prb = hcat([prb[:, 1, s] for s in axes(prb, 3)]...) 
    end

    fi1 = hcat([abs.(prb[:, s] .- 0.5) for s in axes(prb, 2)]...)
    
    efi1 = vec(mean(fi1, dims = 2))
    sbj = eachindex(efi1)
    fi = 4 .* (cumsum(efi1) ./ sbj)

    return fi
end