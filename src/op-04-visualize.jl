"""Function used to visualize simulated operational characteristics.

# Call
`plot(op; kwargs...)`

# Arguments
- `op::DataFrame` -- a data frame, containing simulated operational characteristics; each column represents simulated output for a particular randomization procedure.
- `kwargs` refers to the _kew words_. Here, it is possible to pass _key-value_ pairs supported by a `StatsPlots.heatmap` function.

# Result
- A plot of corresponding operational chacteristics vs. allocation step
"""
function plot(op::DataFrame; kwargs...)
    # getting number of subjects
    nsbj = nrow(op)
    
    # getting number of designs
    ndesign = ncol(op)
    
    # setting Y-axis maximum bound
    yB = ceil(maximum([maximum(row) for row in eachrow(op)]))

    # setting colors
    color_scheme = ColorSchemes.tab10
    ncolors = length(color_scheme)
    colors = ncolors >= ndesign ? hcat([[color_scheme[i]] for i in 1:ndesign]...) : :auto

    # setting shapes
    shapes = [:circle :rect :star4 :diamond :star8 :utriangle :star7 :dtriangle :star5 :rtriangle];
    markers = ncolors >= ndesign ? shapes : :auto

    # making a plot 
    @df op StatsPlots.plot(
        cols();
        size = (800, 600),
        dpi = 300,
        legend = :outerright,
        xlabel = "allocation step",
        ylabel = "value", 
        xlims = (0, nsbj+1),
        ylims = (-0.1, yB),
        xticks = [1; 5:5:nsbj],
        color = colors,
        marker = markers, 
        markercolor = colors,
        markersize = 7,
        markerstrokewidth = 0.5, 
        kwargs...
    )
end


"""Function used to visualize balance-randomness trade-off as a heatmap plot.

# Call
`heatmap(brt; kwargs...)`

# Arguments
- `brt::DataFrame` -- a data frame, containing simulated balance-randomness trade-off; each column represents simulated output for a particular randomization procedure.
- `kwargs` refers to the _kew words_. Here, it is possible to pass _key-value_ pairs supported by a `StatsPlots.heatmap` function.

# Result
- A heatmap plot of the balance-randomness trade-off vs. allocation step
"""
function heatmap(brt::DataFrame; kwargs...)
    brt_transfromed = @pipe brt |> 
        insertcols(_, 1, :sbj => 1:nrow(_)) |> 
        stack(_, Not(:sbj), variable_name = :design) |> 
        unstack(_, :sbj, :value) |>
        sort(_, names(_)[end])

    # getting number of subjects
    nsbj = ncol(brt_transfromed)-1
    
    # getting number of designs
    ndesign = nrow(brt_transfromed)

    # getting design names
    designs = brt_transfromed.design

    # setting colors
    colors = cgrad([:red, :orange, :yellow, :green, :blue, :navy, :purple])

    # making a plot
    @df brt_transfromed[:, Not(:design)] StatsPlots.heatmap(
        cols();
        size = (800, 600),
        dpi = 300,
        legend = :outerright,
        xlabel = "allocation step",
        ylabel = "randomization procedure", 
        xticks = [1; 5:5:nsbj],
        yticks = (1:ndesign, designs), 
        color = colors,
        kwargs...
    )
end


"""Function used to visualize a distribution of the final imbalance as a violin plot.

# Call
`violin(final_imb; kwargs...)`

# Arguments
- `final_imb::DataFrame` -- a data frame, containing simulated final imbalances; each column represents simulated output for a particular randomization procedure.
- `kwargs` refers to the _kew words_. Here, it is possible to pass _key-value_ pairs supported by a `StatsPlots.violin` function.

# Result
- A violin plot of the final imbalances.
"""
function violin(final_imb::DataFrame; kwargs...)
    # defining bounds for the X-axis
    xmin = minimum([minimum(row) for row in eachrow(final_imb)])
    xmax = maximum([maximum(row) for row in eachrow(final_imb)])
    B = maximum(abs.([xmin, xmax]))

    # getting number of designs
    ndesign = ncol(final_imb)

    # getting design names
    designs = brt_transfromed.design

    # setting colors
    color_scheme = ColorSchemes.tab10
    ncolors = length(color_scheme)
    colors = ncolors >= ndesign ? hcat([[color_scheme[i]] for i in 1:ndesign]...) : :auto

    # making a plot
    @df final_imb StatsPlots.violin(
        cols(); 
        permute=(:x, :y),
        size = (800, 600),
        dpi = 300,
        legend = nothing,
        ylabel = "imbalance",
        xlabel = "randomization procedure", 
        yticks = -B:2:B,
        xticks = (1:ndesign, designs), 
        yrotation = 45,
        color = colors,
        kwargs...
    )
end