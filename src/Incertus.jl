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
include("rnd-06-abcd.jl")
include("rnd-07-gbcd.jl")
include("rnd-08-bsd.jl")
include("rnd-09-bcdwit.jl")
include("rnd-10-bud.jl")



include("utils.jl")

export CompleteRandomization
export CRD
export RestrictedRandomization
export CRD, PBD, RAND, TBD, EBCD, ABCD, GBCD, BSD, BCDWIT, BUD

export allocation_prb
export set_label

end # module Incertus
