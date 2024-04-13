"""A type of rectricted randomization, representing _**G**_eneralized _**B**_iased _**C**_oin _**D**_esign (_**GBCD**_).

`GBCD(γ)` command initializes _generalized biased coin design_ with a parameter equal to `γ`, 
targeting `1:1` allocation.

An output of the command is an isntance of GBCD."""
struct GBCD <: RestrictedRandomization 
    target::Vector{Int64}
    parameter::Number

    function GBCD(target::Vector{Int64}, parameter::Number)
        # getting number of treatments
        ntrt = length(target)
        
        @assert parameter >= 0 "The procedures parameter, `γ`, must be non-negative";
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";

        return new(target, parameter)
    end
end
GBCD(parameter::Number) = GBCD([1, 1], parameter)


"""Function calculates allocation probabilities for GBCD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::GBCD`: an object, representing Generalized Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::GBCD, N::Vector{Int64})
    # getting a parameter value of the randomization procedure
    γ = rnd.parameter

    # indicating a current allocation step (current subject's ID)
    j = sum(N) + 1

    # probabilities of treatments' assignments
    P = j == 1 ? 0.5 : N[2]^γ/(N[1]^γ + N[2]^γ) 
    prb = [P, 1-P]

    return prb
end