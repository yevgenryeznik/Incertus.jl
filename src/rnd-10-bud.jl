"""A type of rectricted randomzation, representing _**B**_lock _**U**_rn _**D**_esign (_**BUD**_).

`BUD(λ)` command initializes _block urn desin_ with a parameter equal to ``λ``, 
targeting `1:1` allocation.

`BUD(w, λ)` command initializes _block urn desin_ with a parameter equal to ``λ``, 
targeting allocation specified by `w`.

An output of both commands is an instance of `BUD`.
"""
struct BUD <: RestrictedRandomization
    target::Vector{Int64}
    parameter::Int64

    function BUD(target::Vector{Int64}, parameter::Int64) 
        @assert eltype(target) == Int64 "The procedure isn't implemented for non-integer target allocation!";
        @assert parameter > 0 "The procedure's parameter, `λ`, must be positive!";    
        
        return new(simplify(target), parameter)
    end
    
end
BUD(parameter::Int64) = BUD([1, 1], parameter)

"""Function calculates allocation probabilities for BUD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::BUD`: an object, representing Block Urn Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::BUD, N::Vector{Int64})
    # parameter of the randomization procedure (BUD)
    λ = rnd.parameter

    # target allocation
    w = rnd.target
    W = sum(w)
    
    # current subject's ID
    j = sum(N) + 1

    # the number of minimal balanced sets in previous assignments
    k = minimum([floor(N[i]/w[i]) for i in eachindex(w)])

    # probabilities of treatments' assignments
    prb = [(w[i]*(λ + k) - N[i])/(W*(λ + k) - (j-1)) for i in eachindex(w)]

    return prb 
end