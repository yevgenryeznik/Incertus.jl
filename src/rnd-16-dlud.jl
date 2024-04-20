"""A type of rectricted randomization, representing _**D**_rop-the-_**L**_oser _**U**_rn _**D**_esign (_**DLUD**_).

`DLUD(a)` command initializes _drop-the-loser urn design_ with a parameter equal to `a`,
targeting `1:1` allocation.

`DLUD(w, a)` command initializes _drop-the-loser urn design_ with a parameter equal to `a`,
targeting allocation specified by `w`.

An output of the command is an isntance of DLUD.
"""
struct DLUD <: RestrictedRandomization 
    target::Vector{Int64}
    parameter::Number

    function DLUD(target::Vector{Int64}, parameter::Number)    
        @assert parameter > 0 "The procedure's parameter, `a`, must be larger than 0";
  
        return new(target, parameter)
    end
end
DLUD(parameter::Number) = DLUD([1, 1], parameter)


"""Function simulates a randomization procedure

# Call: 
`simulate(rnd, nsbj, nsim, seed)`

# Arguments
- `rnd::DLUD` : an object, representing a DLUD randomization procedure to be simulated.
- `nsbj::Int64`: number of subjects simulated.
- `nsim::Int64`: number of simulations performed.
- `seed::Int64`: a random seed (for reproducibility); a default is set to 314159.

# Result
- an object, representing simlated randomization, an instance of `SimulatedRandomization`.
"""
function simulate(rnd::DLUD, nsbj::Int64, nsim::Int64, seed::Int64 = 314159)
    # target allocation
    w = rnd.target

    # DLD parameter
    a = rnd.parameter

    # getting number of treatments
    ntrt = length(w)

    # an array to collect treatment assignments
    trt = zeros(Int64, nsbj, ntrt, nsim)

    # an array to collect allocation probabilities
    prb = zeros(Float64, nsbj, ntrt, nsim)

    seed!(seed)
    @showprogress for s in 1:nsim
        # numbers of subjects allocated to treatments (treatment numbers)
        N = zeros(Int64, ntrt)

        # initial urn state
        urn = [0; eachindex(w)]     # an urn state, i.e., balls' types in the urn
        nball = [1; w]              # numbers of balls in an urn
        for j in 1:nsbj
            flag = true
            while flag
                # sampling a ball from the urn
                ball = sample(urn, weights(nball))
                if ball == 0
                    nball[2:end] += a .* w
                else
                    # calculating probability of treatment assignment, 
                    # given the current state of the urn
                    prb[j, :, s] = nball[2:end] ./ sum(nball[2:end])

                    # here, the tretament assignment is made
                    trt[j, ball, s] = 1

                    # removing a ball from urn
                    nball[ball+1] -= 1

                    # updating flag
                    flag = false
                end
            end
        end
    end    
    # setting label
    lbl = label(rnd)

    return SimulatedRandomization(w, lbl, trt, prb)
end