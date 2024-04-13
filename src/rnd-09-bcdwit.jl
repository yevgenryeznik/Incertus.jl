"""A type of rectricted randomization, representing _**B**_iased _**C**_oin _**D**_esign 
_**W**_ith _**I**_mbalance _**T**_olerance (_**BCDWIT**_).

`BCDWIT(p, mti)` command initializes _biased coin design with imbalance tolerance_ 
with _parameters_ equal to `p` and `mti`, targeting `1:1` allocation.

An output of the command is an isntance of BCDWIT.
"""
struct BCDWIT <: RestrictedRandomization
    target::Vector{Int64}
    parameter1::Number # probability of a biased coin
    parameter2::Int64  # MTI parameter

    function BCDWIT(target::Vector{Int64}, parameter1::Number, parameter2::Int64)
        # getting number of treatments
        ntrt = length(target)
        
        @assert parameter2 > 0 "The procedure's parameter, `mti`, must be positive";
        @assert 0.5 <= parameter1 <= 1 "The procedure's parameter, `p`, must be between 0.5 and 1";
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";

        return new(target, parameter1, parameter2)
    end
end
BCDWIT(parameter1::Number, parameter2::Int64) = BCDWIT([1, 1], parameter1, parameter2)


"""Function calculates allocation probabilities for BCDWIT, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::BCDWIT`: an object, representing Biased Coin Design With Imbalance Tolerance.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::BCDWIT, N::Vector{Int64})
    # getting a probability of a biased coin
    p = Float64(rnd.parameter1)

    # getting MTI parameter
    mti = rnd.parameter2  

    # calculating current treatment imbalance
    d = N[1] - N[2]
 
    # probabilities of treatments' assignments
    P = abs(d) < mti ? (d == 0 ? 0.5 : (d < 0 ? p : 1-p)) : (d == mti ? 0 : 1)  
    prb = [P, 1-P]

    return prb
end