"""A type of rectricted randomization, representing _**E**_fron's _**B**_iased _**C**_oin _**D**_esign (_**EBCD**_).

`EBCD(p)` command initializes _Efron's biased coin design_ with a parameter equal to `p`,
targeting `1:1` allocation.

An output of the command is an isntance of EBCD.
"""
struct EBCD <: RestrictedRandomization 
    target::Vector{Int64}
    parameter::Number

    function EBCD(target::Vector{Int64}, parameter::Number)
        # getting number of treatments
        ntrt = length(target)
    
        @assert 0.5 <= parameter <= 1 "The procedure's parameter, `p`, must be between 0.5 and 1";
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";

        return new(target, parameter)
    end
end
EBCD(parameter::Number) = EBCD([1, 1], parameter)


"""Function calculates allocation probabilities for EBCD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::EBCD`: an object, representing Efron's Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::EBCD, N::Vector{Int64})
    # randomization parameter
    p = Float64(rnd.parameter)    

    # imbalance
    d = N[1] - N[2]

    # probabilities of treatments' assignments
    P = d == 0 ? 0.5 : (d < 0 ? p : 1-p)
    prb = [P, 1-P]

    return prb
end