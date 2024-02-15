"""Function used to visualize simulated operational characteristics.

# Call
`visualize(op; kw...)`

# Arguments
- `op::DataFrame` -- a data frame, containing simulated operational characteristics; each column represents simulated output for a particular randomization procedure.

# Key words

`kw` refers to the _kew words_. The following are currently supported and have to be provided:

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


"""Function used to visualize balance-randomness trade-off as a heatmap plot.

# Call
`heatmap(op; kw...)`

# Arguments
- `brt::DataFrame` -- a data frame, containing simulated balance-randomness trade-off; each column represents simulated output for a particular randomization procedure.

# Key words

`kw` refers to the _kew words_. The following are currently supported and have to be provided:

- `xlabel::String` -- a text label for ``X`` axis; default is set ot `"allocation step"`.
- `ylabel::String` -- a text label for ``Y`` axis; default is set ot `"design"`.
- `palette::Vector{Symbol}` -- a vector of colors (symbolic representation) to make a continuous color gradient; default is set to `[:red, :orange, :yellow, :green, :blue, :navy, :purple]`.
- `xticks::Vector{<:Number}` -- ticks for ``X`` axis.
- `size::Tuple{Int64, Int64}` -- an output figure size; default is set to `(800, 600)`.
- `dpi::Int64` -- dpi value; default is set to `300`.

# Result
- A heatmap plot of the balance-randomness trade-off vs. allocation step
"""
function heatmap(
    brt::DataFrame;
    xlabel::String = "allocation step",
    ylabel::String = "design",
    palette::Vector{Symbol} = [:red, :orange, :yellow, :green, :blue, :navy, :purple],
    xticks::Vector{<:Number},
    size::Tuple{Int64, Int64} = (800, 600),
    dpi::Int64 = 300
)
    brt_transfromed = @pipe brt |> 
        insertcols(_, 1, :sbj => 1:nrow(_)) |> 
        stack(_, Not(:sbj), variable_name = :design) |> 
        unstack(_, :sbj, :value) |>
        sort(_, names(_)[end])

    # heatmap plot

    gr(
        size = size, 
        dpi = dpi, 
        legend = :outerright,
        foreground_color_legend = nothing
    )

    @df brt_transfromed[:, Not(:design)]  heatmap(cols(), c=cgrad(palette))
    xlabel!(xlabel)
    ylabel!(ylabel)
    xticks!(xticks)
    yticks!(1:nrow(brt_transfromed), brt_transfromed.design)
end