"""A type of rectricted randomization, representing _**P**_ermuted _**B**_lock _**D**_esign (_**PBD**_).

`PBD(λ)` command initializes a _permuted block design_ with a _block size_ 
equal to `2λ`, targeting `1:1` allocation.

`PBD(w, λ)` command initializes a _permuted block design_  with a parameter `λ`, 
targeting allocation specified by `w`; block size equals to `λ*sum(w)`.

An output of both commands is an instance of `PBD`.
"""
struct PBD <: RestrictedRandomization 
    target::Vector{Int64}
    parameter::Int64

    function PBD(target::Vector{Int64}, parameter::Int64) 
        @assert eltype(target) == Int64 "The procedure isn't implemented for non-integer target allocation";

        return new(simplify(target), parameter)
    end
end
PBD(parameter::Int64) = PBD([1, 1], parameter)


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
    λ = rnd.parameter

    # a block size
    bs = λ*sum(w)

    # indicating a current allocation step (current subject's ID)
    j = sum(N) + 1

    # number of complete blocks among the previous assignments
    k = floor((j-1)/bs)

    # probabilities of treatments' assignments
    prb = [(λ*w[i]*(1 + k) - N[i])/(bs*(1 + k) - (j-1)) for i in eachindex(w)]

    return prb
end