"""A type of rectricted randomzation, representing _**B**_ayesian _**B**_iased _**C**_oin _**D**_esign (_**BBCD**_).

`BBCD(γ, n)` command initializes _Bayesian biased coin design_ with a parameter equal to ``γ``, 
targeting `1:1` allocation in a trial with a _sample size_ equal to `n`.

An output of the command is an isntance of BBCD.
"""
struct BBCD <: RestrictedRandomization
    target::Vector{Int64}
    parameter1::Number # parameter γ
    parameter2::Int64  # sample size


    function BBCD(target::Vector{Int64}, parameter1::Number, parameter2::Int64) 
        # getting number of treatments
        ntrt = length(target)
    
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials!";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation!";
        @assert eltype(target) == Int64 "The procedure isn't implemented for non-integer target allocation!";
        @assert parameter1 > 0 "The procedure's parameter, `γ`, must be positive!";    
        
        return new(simplify(target), parameter1, parameter2)
    end
    
end
BBCD(parameter1::Number, parameter2::Int64) = BBCD([1, 1], parameter1, parameter2)

"""Function calculates allocation probabilities for BBCD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::BBCD`: an object, representing Bayesian Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::BBCD, N::Vector{Int64})
    # parameters of the randomization procedure (BBCD)
    γ = rnd.parameter1
    n = rnd.parameter2

    # current subject's ID
    j = N[1] + N[2]

    # probabilities of treatments' assignments
    P = j+1 == 1 ? 0.5 : (j+1 == 2 ? Int(N[1] == 0) : (n + N[2]/N[1])^(1/γ)/((n + N[2]/N[1])^(1/γ) + (n + N[1]/N[2])^(1/γ)))
    prb = [P, 1-P]

    return prb 
end