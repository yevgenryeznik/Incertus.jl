"""A type of rectricted randomization, representing _**E**_fron's _**B**_iased _**C**_oin _**D**_esign (_**EBCD**_).
"""
struct EBCD <: RestrictedRandomization 
    target::Vector{Int64}
    p::Number

    function EBCD(target::Vector{Int64}, p::Number)
        # getting number of treatments
        ntrt = length(target)
    
        if ntrt > 2 
            throw(NotImplementedForMultiArmTrial())
        elseif !allequal(target) 
            throw(NotImplementedForUnequalAllocation())
        elseif !(0.5 <= p <= 1)
            throw(ParamIsOutOfBounds(:p, 0.5, 1))
        else
            return new([1, 1], p)
        end
    end
end
EBCD(p::Number) = EBCD([1, 1], p)


# function calculates allocation probability for Efron's Biased Coin Design
"""Function calculates allocation probabilities for EBCD, given treatment numbers.

# Arguments
- `rnd::EBCD`: an object, representing Efron's Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::EBCD, N::Vector{Int64})
    # randomization parameter
    p = rnd.p    

    # imbalance
    d = N[2] - N[1]

    # probability of treatment assignment
    return d == 0 ? 0.5 : (d < 0 ? p : 1-p)
end