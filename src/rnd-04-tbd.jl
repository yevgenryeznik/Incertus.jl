"""A type of rectricted randomization, representing _**T**_runcated _**B**_inomial _**D**_esign (_**TBD**_).
"""
struct TBD <: RestrictedRandomization 
    target::Vector{<:Number}
    nsbj::Int64

    function TBD(target::Vector{<:Number}, nsbj::Int64)
        # getting number of treatments
        ntrt = length(target)
        
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";
            
        return new([1, 1], nsbj)
    end
end
TBD(nsbj::Int64) = TBD([1, 1], nsbj)


"""Function calculates allocation probabilities for TBD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::TBD`: an object, representing Truncated Binomial Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::TBD, N::Vector{Int64})
    # sample size
    nsbj = rnd.nsbj

    # allocation probability, given current treatment assignments
    prb = max(N[1], N[2]) < nsbj/2 ? 0.5 : (N[1] < N[2] ? 1 : 0)

    return prb
end