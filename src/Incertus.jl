module Incertus

using DataFrames
using Distributions
using Match
using Pipe
using ProgressMeter
using Random: seed!
using Roots
using StatsPlots: @df, plot, gr, xticks!, yticks!, xlabel!, ylabel!, xlims!, ylims!
using Statistics
 
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
include("rnd-11-eud.jl")
include("rnd-12-bbcd.jl")

include("simulation.jl")
include("op-01-imb.jl")
include("op-02-rand.jl")
include("op-03-brt.jl")
include("op-04-visualize.jl")
include("utils.jl")

export DataFrame, nrow

export CompleteRandomization
export CRD
export RestrictedRandomization
export CRD, PBD, RAND, TBD, EBCD, ABCD, GBCD, BSD, BCDWIT, BUD, EUD, BBCD
export SimulatedRandomization

export simulate
export allocation_prb

export calc_expected_abs_imb
export calc_variance_of_imb
export calc_expected_max_abs_imb
export calc_cummean_loss
export calc_cummean_epcg
export calc_cummean_pda
export calc_fi
export calc_brt
export visualize
export heatmap 

export set_label

end # module Incertus
