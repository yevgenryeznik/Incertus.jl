"""A type of rectricted randomization, representing _**M**_ass _**W**_eighted _**U**_rn _**D**_esign (_**MWUD**_).

`MWUD(α)` command initializes _mass weighted urn design_ with a parameter equal to `α`,
targeting `1:1` allocation.

`MWUD(w, α)` command initializes _mass weighted urn design_ with a parameter equal to `α`,
targeting allocation specified by `w`.

An output of the command is an isntance of MWUD.
"""
struct MWUD <: RestrictedRandomization 
    target::Vector{Int64}
    parameter::Number

    function MWUD(target::Vector{Int64}, parameter::Number)    
        @assert parameter > 0 "The procedure's parameter, `α`, must be larger than 0";
  
        return new(target, parameter)
    end
end
MWUD(parameter::Number) = MWUD([1, 1], parameter)


"""Function calculates allocation probabilities for MWUD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::MWUD`: an object, representing Mass Weighted Urn Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::MWUD, N::Vector{Int64})
    # target allocation
    w = rnd.target
    
    # parameter of the randomization procedure (MWUD)
    α = rnd.parameter
    
    # current subject's ID
    j = sum(N) + 1

    # target allocation proportions
    ρ = w ./ sum(w)
    
    num = (α .* ρ - N + (j-1) .* ρ) .* ((α .* ρ - N + (j-1) .* ρ) .> 0)

    # probabilities of tretament assignments
    return num ./ sum(num)
end

