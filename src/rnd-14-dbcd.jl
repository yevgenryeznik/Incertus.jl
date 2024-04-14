
"""A type of rectricted randomization, representing _**D**_oubly-Adaptive _**B**_iased _**C**_oin _**D**_esign (_**DBCD**_).

`DBCD(γ)` command initializes _doubly-adaptive biased coin design_ with a parameter equal to `γ`,
targeting `1:1` allocation.

`DBCD(w, γ)` command initializes _doubly-adaptive biased coin design_ with a parameter equal to `γ`,
targeting allocation specified by `w`.

An output of the command is an isntance of DBCD.
"""
struct DBCD <: RestrictedRandomization 
    target::Vector{Int64}
    parameter::Number

    function DBCD(target::Vector{Int64}, parameter::Number)
        @assert parameter > 0 "The procedure's parameter, `γ`, must be larger than 0";
  
        return new(target, parameter)
    end
end
DBCD(parameter::Number) = DBCD([1, 1], parameter)


"""Function calculates allocation probabilities for DBCD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::DBCD`: an object, representing Doubly-Adaptive Biased Coin Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::DBCD, N::Vector{Int64})
    # target allocation
    w = rnd.target
    
    # parameter of the randomization procedure (MWUD)
    γ = rnd.parameter

    # target allocation proportions
    ρ = w ./ sum(w)

    if ~all(N .> 0)
        # probabilities of tretament assignments are equal to 
        # target allocation proportions
        return ρ
    else
        # current subject's ID
        j = sum(N) + 1

        num = ρ .* (ρ ./ (N./(j-1))).^γ

        # probabilities of tretament assignments
        return num ./ sum(num)
    end
end