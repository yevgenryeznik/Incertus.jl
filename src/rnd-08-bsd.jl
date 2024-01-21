"""A type of rectricted randomization, representing _**B**_ig _**S**_tick _**D**_esign (_**BSD**_).
"""
struct BSD <: RestrictedRandomization
    target::Vector{<:Number}
    mti::Int64

    function BSD(target::Vector{<:Number}, mti::Int64)
        # getting number of treatments
        ntrt = length(target)
    
        if ntrt > 2 
            throw(NotImplementedForMultiArmTrial())
        elseif !allequal(target) 
            throw(NotImplementedForUnequalAllocation())
        else
            return new([1, 1], mti)
        end
    end
end
BSD(mti::Int64) = BSD([1, 1], mti)


"""Function calculates allocation probabilities for BSD, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::BSD`: an object, representing Big Stick Design.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::BSD, N::Vector{Int64})
    # getting MTI parameter
    mti = rnd.mti

    # calculating current treatment imbalance
    d = N[1] - N[2]
 
    # calculating probability of tretament asssgnment, given imbalance
    prb = abs(d) < mti ? 0.5 : (d == mti ? 0 : 1)    

    return prb
end