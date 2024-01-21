"""A type of rectricted randomization, representing _**A**_djustable _**B**_iased _**C**_oin _**D**_esign (_**ABCD**_).
"""
struct ABCD <: RestrictedRandomization
    target::Vector{Int64}
    a::Number

    function ABCD(target::Vector{Int64}, a::Number)
        # getting number of treatments
        ntrt = length(target)
        
        @assert a >= 0 "The procedures parameter, `a`, must be non-negative";
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";
        
        return new(target, a)
    end
end
ABCD(a::Number) = ABCD([1, 1], a)


"""Function calculates allocation probabilities for ABCD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::ABCD`: an object, representing Adjustable Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::ABCD, N::Vector{Int64})
    # parameter of the randomization procedure
    a = rnd.a

    # current treatment imbalance
    d = N[1] - N[2]

    # allocation function proposed by A. Baldi Antognini and A. Giovagnoli
    p = abs(d)^a/(abs(d)^a + 1)

    # allocation probability, given current imbalance
    prb = abs(d) <= 1 ? 0.5 : (d < -1 ? p : 1-p)    

    return prb
end