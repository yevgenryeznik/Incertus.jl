"""A data structure to deal with _**A**_llocation _**R**_atio _**P**_reserving (_**ARP**_) property.

It has the following fields:

- `label::String`: a label for a randomization procedure that has been simulated.
- `ρ::Vector{Float64}`: a vector of target allocation proportions.
- `expected_prb::Matrix{Float64}`: _expected values_ of ``P_k(j)`` evaluated via simulations.

"""
struct ARP
    label::String
    ρ::Vector{Float64}
    expected_prb::Matrix{Float64}
end


function Base.show(io::IO, arp::ARP)
    lbl = arp.label
    ρstr = join([@sprintf("%4.2f", arp.ρ[i]) for i in eachindex(arp.ρ)], ", ")
    prb = arp.expected_prb

    println("Evaluating ARP property for $lbl")
    println("--------------------------------------------------------")
    println("target allocation proportions: (" * ρstr * ")")
    println("unconditional allocation probabilities:")
    display(prb)

    return nothing
end

"""Function evaluates allocation ratio preserving (ARP) property, i.e., calculates _unconditional allocation probabilities_. A procedure has an ARP property if

# Call
`eval_arp(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, representing simulation output.

# Result
- An instance of  `ARP` type, summarizing the _expected values_ of ``P_k(j)`` via simulations.
"""
function eval_arp(sr::SimulatedRandomization)
    # extracting a label of the simulated randomization 
    label = sr.label

    # extracting target allocation of the simulated randomization 
    target = sr.target
    ρ = target ./ sum(target) # converting to target allocation proportions

    # extracting simulated allocation probabilities
    prb = sr.prb

    # taking the mean
    expected_prb = mean(prb, dims = 3)[:, :, 1]

    return ARP(label, ρ, expected_prb)
end


"""Function evaluates allocation ratio preserving (ARP) property, i.e., calculates _unconditional allocation probabilities_. 

# Call
`eval_arp(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `Vector{ARP}` object -- a vector of instances of `ARP`, each summarizing _expected values_ of ``P_k(j)`` via simulations.
"""
function eval_arp(sr::Vector{SimulatedRandomization})
    return [eval_arp(item) for item in sr]
end