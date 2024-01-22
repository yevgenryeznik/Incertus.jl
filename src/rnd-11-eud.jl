"""A type of rectricted randomzation, representing _**E**_hrenfest _**U**_rn _**D**_esign (_**BUD**_).

A command `EUD(mti)` initializes a randomization procedure.

# Input
- An integer number `mti` (``> 0``) is a parameter of the randomization procedure.

# Output
- An instance of `EUD`.
"""
struct EUD <: RestrictedRandomization
    target::Vector{Int64}
    mti::Int64

    function EUD(target::Vector{Int64}, mti::Int64) 
        # getting number of treatments
        ntrt = length(target)
    
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials!";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation!";
        @assert eltype(target) == Int64 "The procedure isn't implemented for non-integer target allocation!";
        @assert mti > 0 "The procedure's parameter, `mti`, must be positive!";    
        
        return new(simplify(target), mti)
    end
    
end
EUD(mti::Int64) = EUD([1, 1], mti)

"""Function calculates allocation probabilities for EUD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::EUD`: an object, representing Ehrenfest Urn Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::EUD, N::Vector{Int64})
    # parameter of the randomization procedure (EUD)
    mti = rnd.mti

    # current imbalance
    d = N[1] - N[2]

    # probability of treatment assignment
    return 0.5*(1 - d/mti) 
end