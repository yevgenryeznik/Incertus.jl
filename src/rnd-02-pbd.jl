"""A type of rectricted randomization, representing _**P**_ermuted _**B**_lock _**D**_esign (_**PBD**_).

`PBD(λ)` command initializes a _permuted block design_ with a _block size_ 
equal to `2λ`, targeting `1:1` allocation.

`PBD(w, λ)` command initializes a _permuted block design_  with a parameter `λ`, 
targeting allocation specified by `w`; block size equals to `λ*sum(w)`.

An output of both commands is an instance of `PBD`.
"""
struct PBD <: RestrictedRandomization 
    target::Vector{Int64}
    param::Int64

    function PBD(target::Vector{Int64}, param::Int64) 
        @assert eltype(target) == Int64 "The procedure isn't implemented for non-integer target allocation";

        return new(simplify(target), param)
    end
end
PBD(param::Int64) = PBD([1, 1], param)


"""Function calculates allocation probabilities for PBD, given treatment numbers.
# Call
- allocation_prb(rnd, N)

# Arguments
- `rnd::PBD`: an object, representing Permuted Block Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::PBD, N::Vector{Int64})
    # target allocation
    w = rnd.target

    # parameter of the randomization procedure
    λ = rnd.param

    # a block size
    bs = λ*sum(w)

    # indicating a current allocation step (current subject's ID)
    j = sum(N) + 1

    # number of complete blocks among the previous assignments
    k = floor((j-1)/bs)

    # number of treatments
    ntrt = length(w)

    # if it is a randomization, targeting 1:1 allocation
    if ntrt == 2 && allequal(w)
        # calculating probability of treatment assignment
        prb = (0.5*bs*(1 + k) - N[1])/(bs*(1 + k)-(j-1))

    # otherwise, it is a randomization, targeting allocation w 
    # (unequal in a two-arm trial or equal/unequal in a multi-arm trial)   
    else
        # vector of probabilities of treatment assignments
        prb = [(λ*w[i]*(1 + k) - N[i])/(bs*(1 + k) - (j-1)) for i in eachindex(w)]
    end  

    return prb
end