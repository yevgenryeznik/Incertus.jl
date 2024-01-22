"""A type of rectricted randomzation, representing _**B**_lock _**U**_rn _**D**_esign (_**BUD**_).

A command `BUD(w, λ)` initializes a randomization procedure.

# Input
- An input `w` represents the target allocation ratio and has a type `Vector{Int64}`.
- An integer number `λ` (``> 0``) is a parameter of the randomization procedure.

# Output
- An instance of `BUD`.

A command `BUD(λ)` initializes a randomization procedure, targeting 1:1 allocation.
"""
struct BUD <: RestrictedRandomization
    target::Vector{Int64}
    λ::Int64

    function BUD(target::Vector{Int64}, λ::Int64) 
        @assert eltype(target) == Int64 "The procedure isn't implemented for non-integer target allocation!";
        @assert λ > 0 "The procedure's parameter, `λ`, must be positive!";    
        
        return new(simplify(target), λ)
    end
    
end
BUD(λ::Int64) = new([1, 1], λ)

"""Function calculates allocation probabilities for BUD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::BUD`: an object, representing Block Urn Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::BUD, N::Vector{Int64})
    # parameter of the randomization procedure (BUD)
    λ = rnd.param

    # target allocation
    w = rnd.target
    W = sum(w)
    
    # number of treatments
    ntrt = length(w)

    # current subject's ID
    j = sum(N) + 1

    if ntrt == 2
        prb = (λ + minimum(N) - N[1])/(2*(λ + minimum(N)) - (j-1))
    else
        # the number of minimal balanced sets in previous assignments
        k = minimum([floor(N[i]/w[i]) for i in eachindex(w)])
    
        # probabilities of tretament assignments
        [(w[i]*(λ + k) - N[i])/(W*(λ + k) - (j-1)) for i in eachindex(w)]
    end

    return prb 
end