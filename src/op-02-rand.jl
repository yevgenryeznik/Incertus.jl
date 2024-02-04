# calculating correct guess probability under the convergence guessing strategy
function calc_guess_cgs(δ::Vector{Int64})
    imb = [0; 2 .* cumsum(δ[1:end-1]) - eachindex(δ[1:end-1])]
    guess = [0.5*(1-sign(item)) for item in imb]

    return(guess)
end


# calculating correct guess probability under the maximum probability guessing strategy
function calc_guess_mpgs(ϕ::Vector{Float64})
    guess = [0.5*(1-sign(item-0.5)) for item in ϕ]

    return(guess)
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
    ntrt = maximum(ntrt)
    if ntrt == 2
        trt = 2 .- trt
    end

    pcg = gs == "C" ? 
        hcat([calc_guess_cgs(trt[:, s]) for s in axes(trt, 2)]...) :
        hcat([calc_guess_mpgs(prb[:, s]) for s in axes(trt, 2)]...)
    
    sbj = eachindex(pcg)
    expected_pcg = vec(mean(abs_imb, dims = 2))
    cummean_of_epcg = cumsum(expected_pcg) ./ sbj

    return cummean_of_epcg
end