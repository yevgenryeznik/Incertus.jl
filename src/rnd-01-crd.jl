"""A type of randomization, representing _**C**_ompletely _**R**_andomized _**D**_esign (_**CRD**_).

`CRD()` command initializes a complete randomization procedure, targeting `1:1` allocation.

`CRD(w)` command initializes a complete randomization procedure, targeting allocation specified by `w`.

An output of both commands is an instance of `CRD`.
"""
struct CRD <: CompleteRandomization
    target::Vector{<:Number}
    CRD(target::Vector{<:Number}) = new(simplify(target))
end
CRD() = CRD([1, 1])


"""Function calculates allocation probabilities for CRD.
# Call
- `allocation_prb(rnd)`

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