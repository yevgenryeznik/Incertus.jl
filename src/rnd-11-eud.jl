"""A type of rectricted randomzation, representing _**E**_hrenfest _**U**_rn _**D**_esign (_**BUD**_).

`EUD(mti)` command initializes _Ehrenfest's urn design_ with a parameter equal to ``mti``, 
targeting `1:1` allocation.

An output of the command is an isntance of EUD.
"""
struct EUD <: RestrictedRandomization
    target::Vector{Int64}
    parameter::Int64 # MTI parameter

    function EUD(target::Vector{Int64}, parameter::Int64) 
        # getting number of treatments
        ntrt = length(target)
    
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials!";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation!";
        @assert eltype(target) == Int64 "The procedure isn't implemented for non-integer target allocation!";
        @assert parameter > 0 "The procedure's parameter, `mti`, must be positive!";    
        
        return new(simplify(target), parameter)
    end
    
end
EUD(parameter::Int64) = EUD([1, 1], parameter)

"""Function calculates allocation probabilities for EUD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::EUD`: an object, representing Ehrenfest Urn Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::EUD, N::Vector{Int64})
    # parameter of the randomization procedure (EUD)
    mti = rnd.parameter

    # current imbalance
    d = N[1] - N[2]

    # probabilities of treatments' assignments
    P = 0.5*(1 - d/mti)
    prb = [P, 1-P]

    return prb 
end