"""A type of randomization, representing _**C**_ompletely _**R**_andomized _**D**_esign (_**CRD**_).
"""
struct CRD <: CompleteRandomization
    target::Vector{<:Number}
    CRD(target::Vector{<:Number}) = new(simplify(target))
end
CRD() = CRD([1, 1])


# function calculates allocation probability for Complete Randomization
"""Function calculates allocation probabilities for CRD.

# Arguments
- `rnd::CRD`: an object, representing Complete Randomization.
"""
function allocation_prb(rnd::CRD)
    # target alloaction
    w = rnd.target    

    # probability(ies) of treatment assignment
    if allequal(w) && length(w) == 2
        prb = 0.5
    else
        prb = w/sum(w)
    end

    return prb
end