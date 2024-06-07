"""A type of rectricted randomzation, representing _**B**_ayesian _**B**_iased _**C**_oin _**D**_esign (_**BBCD**_).

`BBCD(γ)` command initializes _Bayesian biased coin design_ with a parameter equal to ``γ``, 
targeting `1:1` allocation.

An output of the command is an isntance of BBCD.
"""
struct BBCD <: RestrictedRandomization
    target::Vector{Int64}
    parameter::Number # parameter γ

    function BBCD(target::Vector{Int64}, parameter::Number) 
        # getting number of treatments
        ntrt = length(target)
    
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials!";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation!";
        @assert eltype(target) == Int64 "The procedure isn't implemented for non-integer target allocation!";
        @assert parameter > 0 "The procedure's parameter, `γ`, must be positive!";    
        
        return new(simplify(target), parameter)
    end
    
end
BBCD(parameter::Number) = BBCD([1, 1], parameter)

"""Function calculates allocation probabilities for BBCD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::BBCD`: an object, representing Bayesian Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::BBCD, N::Vector{Int64})
    # parameters of the randomization procedure (BBCD)
    γ = rnd.parameter
    
    # current subject's ID
    j = N[1] + N[2]

    # probabilities of treatments' assignments
    P = 0.5
    if (j+1 == 2)
        P = Int(N[1] == 0)
    end
    if (j + 1 > 2)
        r = N[2]/N[1]
        P = (1 + r/j)^(1/γ)/((1 + r/j)^(1/γ) + (1 + (1/r)/j)^(1/γ))
    end
    prb = [P, 1-P]

    return prb 
end