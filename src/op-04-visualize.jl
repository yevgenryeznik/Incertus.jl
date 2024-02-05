"""Function used to visualize simulated operational characteristics.

# Call
`visualize(op, xlabel, ylabel, colors, shapes, marker_size, xticks, yticks, xlims, ylims, size, dpi)`

# Arguments
- `op::DataFrame` -- a data frame, containing simulated operational characteristics; each column represents simulated output for a particular randomization procedure.
- `xlabel::String` -- a text label for ``X`` axis.
- `ylabel::String` -- a text label for ``Y`` axis.
- `colors::Matrix{Symbol}` -- a matrix (row) of colors (symbolic representation) to distinguish outputs for different randomization procedure.
- `shapes::Matrix{Symbol}` -- a matrix (row) of markers (symbolic representation) to distinguish outputs for different randomization procedure.
- `marker_size::Int64` -- marker size.
- `xticks::Vector{<:Number}` -- ticks for ``X`` axis.
- `yticks::Vector{<:Number}` -- ticks for ``Y`` axis.
- `xlims::Tuple{<:Number, <:Number}` -- a tuple to set limits of ``X`` axis; default is set to `(0, nrow(op)+1)`.
- `ylims::Tuple{<:Number, <:Number}` -- a tuple to set limits of ``Y`` axis; default is set to `(0, ceil(maximum([maximum(row) for row in eachrow(op)])))`.
- `size::Tuple{Int64, Int64}` -- an output figure size; default is set to `(800, 600)`.
- `dpi::Int64` -- dpi value; default is set to `300`.

# Result
- A plot of corresponding operational chacteristics vs. allocation step
"""
function visualize(
    op::DataFrame;
    xlabel::String,
    ylabel::String,
    colors::Matrix{Symbol},
    shapes::Matrix{Symbol},
    marker_size::Int64,
    xticks::Vector{<:Number},
    yticks::Vector{<:Number},
    xlims::Tuple{<:Number, <:Number} = (0, nrow(op)+1),
    ylims::Tuple{<:Number, <:Number} = (0, ceil(maximum([maximum(row) for row in eachrow(op)]))),
    size::Tuple{Int64, Int64} = (800, 600),
    dpi::Int64 = 300
    )

    gr(
        size = size, 
        dpi = dpi, 
        legend = :outerright,
        foreground_color_legend = nothing
    )
    
    @df op plot(cols(), marker = shapes, ms = marker_size, color = colors)

    xticks!(xticks)
    yticks!(yticks)
    xlabel!(xlabel)
    ylabel!(ylabel)
    xlims!(xlims)
    ylims!(ylims)
end