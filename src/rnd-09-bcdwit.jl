"""A type of rectricted randomization, representing _**B**_iased _**C**_oin _**D**_esign 
_**W**_ith _**I**_mbalance _**T**_olerance (_**BCDWIT**_).
"""
struct BCDWIT <: RestrictedRandomization
    target::Vector{Int64}
    p::Number
    mti::Int64

    function BCDWIT(target::Vector{Int64}, p::Number, mti::Int64)
        # getting number of treatments
        ntrt = length(target)
        
        @assert mti > 0 "The procedures parameter, `mti`, must be positive";
        @assert 0.5 <= p <= 1 "The procedures parameter, `p`, must be between 0.5 and 1";
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";


        return new([1, 1], p, mti)
    end
end
BCDWIT(p::Number, mti::Int64) = BCDWIT([1, 1], p, mti)


"""Function calculates allocation probabilities for BCDWIT, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::BCDWIT`: an object, representing Biased Coin Design With Imbalance Tolerance.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::BCDWIT, N::Vector{Int64})
    # getting a probability of a biased coin
    p = Float64(rnd.p)

    # getting MTI parameter
    mti = rnd.mti  

    # calculating current treatment imbalance
    d = N[1] - N[2]
 
    # calculating probability of tretament asssgnment, given imbalance
    prb = abs(d) < mti ? (d == 0 ? 0.5 : (d < 0 ? p : 1-p)) : (d == mti ? 0 : 1)  

    return prb
end