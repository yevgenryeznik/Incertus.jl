"""A type of rectricted randomization, representing _**P**_ermuted _**B**_lock _**D**_esign (_**PBD**_).
"""
struct PBD <: RestrictedRandomization 
    target::Vector{Int64}
    param::Int64

    function PBD(target::Vector{Int64}, param::Int64) 
        return new(simplify(target), param)
    end
end
PBD(param::Int64) = PBD([1, 1], param)


# function calculates allocation probability for Permuted Block Design.
"""Function calculates allocation probabilities for PBD, given treatment numbers.

# Arguments
- `rnd::PBD`: an object, representing Permuted Block Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::PBD, N::Vector{Int64})
    # target allocation
    w = rnd.target

    # number of treatments
    ntrt = length(w)

    # indicating a current allocation step (current subject's ID)
    j = sum(N) + 1

    # if it is a randomization, targeting 1:1 allocation
    if ntrt == 2 && allequal(w)
        # a block size
        bs = 2*rnd.param

        # indicating a current subject's position in a block:
        # 1 <= k <= bs
        k = (j % bs == 0)*bs + (j % bs)

        # calculating the number of subjects in a block 
        # allocated to the 1st treatment arm (n1)
        n1 = k == 1 ? 0 : N[1] - ((j - 1) ÷ bs)*(bs/2) 

        # calculating probability of treatment assignment
        prb = (0.5*bs-n1)/(bs-k+1)
    
    # otherwise, it is a randomization, targeting allocation w 
    # (unequal in a two-arm trial or equal/unequal in a multi-arm trial)   
    else
        # getting procedure's parameter
        λ = rnd.param

        # calculating block size
        bs = λ*sum(w)

        # number of complete blocks among the previous assignments
        k = floor((j-1)/bs)
    
        # vector of probabilities of treatment assignments
        prb = [(λ*w[i]*(1 + k) - N[i])/(bs*(1 + k) - (j-1)) for i in eachindex(w)]
    end  

    return prb
end