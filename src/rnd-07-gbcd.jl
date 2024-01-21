"""A type of rectricted randomization, representing _**G**_eneralized _**B**_iased _**C**_oin _**D**_esign (_**GBCD**_).
"""
struct GBCD <: RestrictedRandomization 
    target::Vector{Int64}
    γ::Number

    function GBCD(target::Vector{Int64}, γ::Number)
        # getting number of treatments
        ntrt = length(target)
    
        if ntrt > 2 
            throw(NotImplementedForMultiArmTrial())
        elseif !allequal(target) 
            throw(NotImplementedForUnequalAllocation())
        else
            return new([1, 1], γ)
        end
    end
end
GBCD(γ::Number) = GBCD([1, 1], γ)


"""Function calculates allocation probabilities for GBCD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::GBCD`: an object, representing Generalized Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function calculate_prb(rnd::GBCD, N::Vector{Int64})
    # getting a parameter value of the randomization procedure
    γ = rnd.γ

    # indicating a current allocation step (current subject's ID)
    j = sum(N) + 1

    # calculating probability of treatment asssgnment
    prb = j == 1 ? 0.5 : N2^γ/(N1^γ + N2^γ) 

    return prb
end