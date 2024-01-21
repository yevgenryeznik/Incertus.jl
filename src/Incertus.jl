module Incertus

using Distributions
using Match
using Plots
using ProgressMeter
using Random
using Roots
 
abstract type CompleteRandomization end
abstract type RestrictedRandomization end

include("rnd-01-crd.jl")
include("rnd-02-pbd.jl")
include("rnd-03-rand.jl")
include("rnd-04-tbd.jl")
include("rnd-05-ebcd.jl")
include("abcd.jl")



include("utils.jl")

export CompleteRandomization
export CRD
export RestrictedRandomization
export ABCD, CRD, EBCD, PBD, RAND, TBD

export allocation_prb
export set_label

end # module Incertus
