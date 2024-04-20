"""A type of rectricted randomization, representing _**T**_runcated _**M**_ultiinomial _**D**_esign (_**TMD**_).

`TMD(n)` command initializes a _truncated multinomial design_, targeting `1:1` allocation
in a trial with a _sample size_ equal to `n`. (In this case, it is similar to `TBD`).

`TMD(w, n)` command initializes a _truncated multinomial design_, targeting allocation specified by `w`
in a trial with a _sample size_ equal to `n`.

An output of the command is an isntance of TMD.
"""
struct TMD <: RestrictedRandomization 
    target::Vector{<:Number}
    parameter::Int64

    function TMD(target::Vector{<:Number}, parameter::Int64)            
        return new(target, parameter)
    end
end
TMD(parameter::Int64) = TMD([1, 1], parameter)


"""Function calculates allocation probabilities for TMD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::TMD`: an object, representing Truncated Multinomial Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::TMD, N::Vector{Int64})
    # target allocation
    w = rnd.target

    # target allocation propotions
    ρ = w ./ sum(w)

    # number of treatments
    ntrt = length(w)

    # sample size
    nsbj = rnd.parameter

    # probabilities of treatments' assignments
    prb = zeros(Float64, ntrt)
    idx = N - nsbj .* ρ .< 0
    prb[idx] = w[idx] ./ sum(w[idx])

    return prb
end