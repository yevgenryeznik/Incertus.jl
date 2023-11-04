"""A type of rectricted randomization, representing _**Rand**_om Allocation Rule (_**Rand**_).
"""
struct RAND <: RestrictedRandomization
    target::Vector{Int64}
    nsbj::Int64

    RAND(target::Vector{Int64}, nsbj::Int64) = new(simplify(target), nsbj)
end
RAND(nsbj::Int64) = RAND([1, 1], nsbj)