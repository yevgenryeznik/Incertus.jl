module Incertus

using CategoricalArrays
using ColorSchemes
using DataFrames
using Distributions
using Latexify
using Match
using NonlinearSolve
using Pipe
using Plots
using Printf
using ProgressMeter
using Random: seed!
using Roots
using StatsBase
using StatsPlots
using Statistics

abstract type Randomization end
abstract type CompleteRandomization <: Randomization end
abstract type RestrictedRandomization <: Randomization end

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
include("rnd-11-eud.jl")
include("rnd-12-bbcd.jl")
include("rnd-13-mwud.jl")
include("rnd-14-dbcd.jl")
include("rnd-15-maxent.jl")
include("rnd-16-dlud.jl")
include("rnd-17-tmd.jl")

include("simulation.jl")
include("op-01-imb.jl")
include("op-02-rand.jl")
include("op-03-brt.jl")
include("op-04-arp.jl")
include("op-05-visualize.jl")
include("utils.jl")

export DataFrame

export Randomization
export CompleteRandomization
export CRD
export RestrictedRandomization
export CRD, PBD, RAND, TBD, EBCD, ABCD, GBCD, BSD, BCDWIT, BUD, EUD, BBCD, MWUD, DBCD, MaxEnt, DLUD, TMD
export SimulatedRandomization
export ARP
export simulate
export allocation_prb

export calc_final_imb
export calc_expected_abs_imb
export calc_variance_of_imb
export calc_expected_max_abs_imb
export calc_cummean_loss
export calc_cummean_epcg
export calc_cummean_pda
export calc_fi
export calc_brt
export eval_arp

export plot
export heatmap
export violin

export label

end # module Incertus
