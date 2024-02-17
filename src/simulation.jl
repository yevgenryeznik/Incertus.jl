"""A type, representing an output of a simulated randomization procedure.

A command `SimulatedRandomization(label, trt, prb)` initializes an instance of `SimulatedRandomization`:

- `label` is a string that describes the simulated randomization procedure; 
- `trt` is a matrix of size ``n\\times S``, representing treatment assignments;
- `prb` is an array of size ``n\\times K\\times S``, representing allocation probabilities,

where 

- ``n`` is the number of subjects simulated;
- ``K`` is the number of treatments simulated (equal to the length of the target allocation vector);
- ``S`` is the number of simulations performed.
"""
struct SimulatedRandomization
    label::String
    trt::Matrix{Int64}
    prb::Array{Float64} 
end

# function used to display simulated randomization output.
function Base.show(io::IO, sr::SimulatedRandomization)
    # getting the description of the simulated randomization procedure
    label = sr.label

    # extracting simulated treatment assignments
    trt = sr.trt

    # extracting simulated allocaton probabilities 
    prb = sr.prb

    println(io, "An output of the randomization procedure ($(label)) simulation (`trt`, `prb`)")
    println(io, "--------------------------------------------------------")
    println(io, "Simulated treatment assignments (`trt`):")
    display(trt)
    println(io, "Simulated allocation probabilities (`prb`):")
    display(prb) 
end



"""Function simulates a randomization procedure

# Call: 
`simulate(rnd, nsbj, nsim, seed)`

# Arguments
- `rnd::T where T <: Randomization`: an object, representing a randomization procedure to be simulated.
- `nsbj::Int64`: number of subjects simulated.
- `nsim::Int64`: number of simulations performed.
- `seed::Int64`: a random seed (for reproducibility); a default is set to 314159.

# Result
- an object, representing simlated randomization, an instance of `SimulatedRandomization`.
"""
function simulate(rnd::T, nsbj::Int64, nsim::Int64, seed::Int64 = 314159) where T <: Randomization
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
    # setting label
    label = set_label(rnd)

    return SimulatedRandomization(label, trt, prb)
end


"""Function simulates a set of randomization procedures

# Call: 
`simulate(rnd, nsbj, nsim, seed)`

# Arguments
- `rnd::Vector{<:Randomization}`: a vector of instances of `<:Randomization` type; each instance represents a randomization procedure to be simulated.
- `nsbj::Int64`: number of subjects simulated.
- `nsim::Int64`: number of simulations performed.
- `seed::Int64`: a random seed (for reproducibility); a default is set to 314159.

# Result
- a vector of instances of `SimulatedRandomization` type; each instance represents a simulated randomization procedure.
"""
function simulate(rnd::Vector{<:Randomization}, nsbj::Int64, nsim::Int64, seed::Int64 = 314159)
    return [simulate(item, nsbj, nsim, seed) for item in rnd]
end