"""A type of rectricted randomization, representing _**Rand**_om Allocation Rule (_**Rand**_).

`RAND(n)` command initializes a _random allocation rule_ , targeting `1:1` allocation 
in a trial with a _sample size_ equal to `n`.

`RAND(w, n)` command initializes a _random allocation rule_ , targeting allocation 
specified by `w` in a trial with a _sample size_ equal to `n`.

An output of both commands is an instance of `RAND`.
"""
struct RAND <: RestrictedRandomization
    target::Vector{Int64}
    parameter::Int64

    RAND(target::Vector{Int64}, parameter::Int64) = new(simplify(target), parameter)
end
RAND(parameter::Int64) = RAND([1, 1], parameter)


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
    nsbj = rnd.parameter

    # indicating a current allocation step (current subject's ID)
    j = sum(N) + 1
    
    # probabilities of treatments' assignments
    W = sum(w)
    prb = [(w[i]/W*nsbj - N[i])/(nsbj - (j-1)) for i in eachindex(w)]

    return prb
end