"""A type of rectricted randomization, representing _**B**_ig _**S**_tick _**D**_esign (_**BSD**_).

`BSD(mti)` command initializes _big stick design_ with a parameter equal to `mti`, 
targeting `1:1` allocation.

An output of the command is an isntance of BSD.
"""
struct BSD <: RestrictedRandomization
    target::Vector{<:Number}
    parameter::Int64

    function BSD(target::Vector{<:Number}, parameter::Int64)
        # getting number of treatments
        ntrt = length(target)
        
        @assert parameter >= 0 "The procedures parameter, `mti`, must be non-negative";
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";

        return new(target, parameter)
    end
end
BSD(parameter::Int64) = BSD([1, 1], parameter)


"""Function calculates allocation probabilities for BSD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::BSD`: an object, representing Big Stick Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::BSD, N::Vector{Int64})
    # getting MTI parameter
    mti = rnd.parameter

    # calculating current treatment imbalance
    d = N[1] - N[2]
 
    # probabilities of treatments' assignments
    P = abs(d) < mti ? 0.5 : (d == mti ? 0 : 1)    
    prb = [P, 1-P]
    
    return prb
end