"""A type, representing an output of a simulated randomization procedure.

A command `SimulatedRandomization(trt, prb)` initializes an instance of `SimulatedRandomization`:

- `trt` is a matrix of size ``n\\times S``, representing treatment assignments;
- `prb` is an array of size ``n\\times K\\times S``, representing allocation probabilities,

where 

- ``n`` is the number of subjects simulated;
- ``K`` is the number of treatments simulated (equal to the length of the target allocation vector);
- ``S`` is the number of simulations performed.
"""
struct SimulatedRandomization
    trt::Matrix{Int64}
    prb::Array{Float64} 
end

# function used to display simulated randomization output.
function Base.show(io::IO, sr::SimulatedRandomization)
    # extracting simulated treatment assignments
    trt = sr.trt

    # extracting simulated allocaton probabilities 
    prb = sr.prb

    println(io, "An output of the randomization simulation (`trt`, `prb`)")
    println(io, "--------------------------------------------------------")
    println(io, "A simulated treatment assignments (`trt`):")
    display(trt)
    println(io, "A simulated allocation probabilities (`prb`):")
    display(prb) 
end



"""Function simulates randomization procedures
`simulate(rnd, nsbj, nsim, seed)`

# Arguments
- `rnd::Union{CompleteRandomization, RestrictedRandomization}`: an object, representing 
a randomization procedure to be simulated, an instance of `CompleteRandomization` or 
`RestrictedRandomization`.
- `nsbj::Int64`: number of subjects simulated.
- `nsim::Int64`: number of simulations performed.
- `seed::Int64`: a random seed (for reproducibility); a default is set to 314159.

# Result
- an object, representing simlated randomization, an instance of `Simulatedrandomization`.
"""
function simulate(rnd::Union{CompleteRandomization, RestrictedRandomization}, nsbj::Int64, nsim::Int64, seed::Int64 = 314159)
    # setting random seed (for reproducibility)
    seed!(seed)

    # target allocation ratio 
    w = rnd.target

    # number of treatments 
    ntrt = length(w)
    
    # 2D-array of treatment assignments
    trt = zeros(Int64, nsbj, nsim)

    # 3D-array of probabilities of treatment assignments
    prb = zeros(Float64, nsbj, ntrt, nsim)
    @showprogress for s in 1:nsim
        # initializing treatment numbers
        N = zeros(Int64, ntrt)
        for j in 1:nsbj
            # probability of treatment assignemnt
            p = typeof(rnd) == CRD ? allocation_prb(rnd) : allocation_prb(rnd, N)
            prb[j, :, s] = ntrt == 2 ? [p, 1-p] : p

            # making treatment assignment
            k = rand(Categorical(prb[j, :, s]))
            trt[j, s] = k

            # updating treatment numbers
            N[k] += 1
        end
    end

    return SimulatedRandomization(trt, prb)
end