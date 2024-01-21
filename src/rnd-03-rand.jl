"""A type of rectricted randomization, representing _**Rand**_om Allocation Rule (_**Rand**_).
"""
struct RAND <: RestrictedRandomization
    target::Vector{Int64}
    nsbj::Int64

    RAND(target::Vector{Int64}, nsbj::Int64) = new(simplify(target), nsbj)
end
RAND(nsbj::Int64) = RAND([1, 1], nsbj)


"""Function calculates allocation probabilities for RAND, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::RAND`: an object, representing Random Allocation Rule.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::RAND, N::Vector{Int64})
    # target allocation
    w = rnd.target
        
    # getting a block size (for Rand, it is equal to the total sample size)
    nsbj = rnd.nsbj

    # indicating a current allocation step (current subject's ID)
    j = sum(N) + 1

    # number of treatments
    ntrt = length(w)
    
    # if it is a randomization, targeting 1:1 allocation
    if ntrt == 2 && allequal(w)
        # calculating probability of treatment assignment
        prb = (0.5*nsbj-N[1])/(nsbj-(j-1))

    # otherwise, it is a randomization, targeting allocation w 
    # (unequal in a two-arm trial or equal/unequal in a multi-arm trial)   
    else
        W = sum(w)

        # vector of probabilities of treatment assignments
        prb = [(w[i]/W*nsbj - N[i])/(nsbj - (j-1)) for i in eachindex(w)]
    end  

    return prb
end