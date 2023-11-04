module Incertus

using Distributions
using Match
using Plots
using ProgressMeter
using Random
using Roots
 
abstract type CompleteRandomization end
abstract type RestrictedRandomization end

include("abcd.jl")
include("crd.jl")
include("ebcd.jl")
include("pbd.jl")
include("tbd.jl")

include("utils.jl")

export CompleteRandomization
export CRD
export RestrictedRandomization
export ABCD, CRD, EBCD, PBD, TBD

export allocation_prb
export set_label

end # module Incertus
