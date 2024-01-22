"""A type of rectricted randomization, representing _**E**_fron's _**B**_iased _**C**_oin _**D**_esign (_**EBCD**_).

`EBCD(p)` command initializes _Efron's biased coin design_ randomization procedure, 
targeting `1:1` allocation in a trial with a _parameter_ equal to ``p``.

An output of the command is an isntance of EBCD.
"""
struct EBCD <: RestrictedRandomization 
    target::Vector{Int64}
    p::Number

    function EBCD(target::Vector{Int64}, p::Number)
        # getting number of treatments
        ntrt = length(target)
    
        @assert 0.5 <= p <= 1 "The procedures parameter, `p`, must be between 0.5 and 1";
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";

        return new(target, p)
    end
end
EBCD(p::Number) = EBCD([1, 1], p)


"""Function calculates allocation probabilities for EBCD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::EBCD`: an object, representing Efron's Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::EBCD, N::Vector{Int64})
    # randomization parameter
    p = Float64(rnd.p)    

    # imbalance
    d = N[2] - N[1]

    # probability of treatment assignment
    return d == 0 ? 0.5 : (d < 0 ? p : 1-p)
end