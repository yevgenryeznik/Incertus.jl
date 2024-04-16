"""Function used to visualize simulated operational characteristics.

# Call
`plot(op; kwargs...)`

# Arguments
- `op::DataFrame` -- a data frame, containing simulated operational characteristics; each column represents simulated output for a particular randomization procedure.
- `kwargs` refers to the _kew words_. Here, it is possible to pass _key-value_ pairs supported by a `StatsPlots.plot` function.

# Result
- A plot of corresponding operational chacteristics vs. allocation step
"""
function plot(op::DataFrame; kwargs...)
    # getting number of subjects
    nsbj = nrow(op)
    
    # getting number of designs
    ndesign = ncol(op)
    
    # setting Y-axis maximum bound
    yB = maximum([maximum(row) for row in eachrow(op)])

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
        legend_foreground_color = nothing,
        xlabel = "allocation step",
        xlims = (0, nsbj+1),
        ylims = (-yB/50, yB*51/50),
        xticks = [1; 5:5:nsbj],
        color = colors,
        marker = markers, 
        markercolor = colors,
        markersize = 7,
        markerstrokewidth = 0.5, 
        kwargs...
    )
end


"""Function used to visualize simulated operational characteristics.

# Call
`plot(op; kwargs...)`

# Arguments
- `op::ARP` -- an instance of `ARP`, containing unconditional allocation probabilities.
- `kwargs` refers to the _kew words_. Here, it is possible to pass _key-value_ pairs supported by a `StatsPlots.plot` function.

# Result
- A plot of corresponding operational chacteristics vs. allocation step
"""
function plot(op::ARP; kwargs...)
    # design label
    lbl = op.label

    # target allocation proportions
    ρ = op.ρ
    
    # unconditional allocation probabilities
    prb = op.expected_prb

    # setting colors
    color_scheme = ColorSchemes.tab10
    ncolors = length(color_scheme)
    colors = ncolors >= size(prb, 2) ? hcat([[color_scheme[i]] for i in axes(prb, 2)]...) : :auto

    # setting shapes
    shapes = [:circle :rect :star4 :diamond :star8 :utriangle :star7 :dtriangle :star5 :rtriangle];
    markers = ncolors >= size(prb, 2) ? shapes : :auto

    # transforming `prb` matrix into a data frame
    prb_df = DataFrame(prb, [latexify("π_$i") for i in axes(prb, 2)])

    # making a plot 
    hline(ρ, ls = :dash, lw = 2.0, lc = :black, label = "target allocation")
    @df prb_df StatsPlots.plot!(
        cols();
        size = (800, 600),
        dpi = 300,
        legend = :topright,
        legend_foreground_color = nothing,
        xlabel = "allocation step",
        ylabel = "unconditional allocation probability",
        title = lbl, 
        xlims = (0, nsbj+1),
        ylims = (0.0, 1.0),
        xticks = :auto,
        yticks = 0:0.1:1,
        color = colors,
        marker = markers, 
        markercolor = colors,
        markersize = 3,
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
    # getting number of simulations
    nsim = nrow(final_imb)

    # getting number of designs
    ndesign = ncol(final_imb)

    # getting design names
    designs = names(final_imb)

    # setting colors
    color_scheme = ColorSchemes.tab10
    ncolors = length(color_scheme)
    colors = [ncolors >= ndesign ? color_scheme[i] : :auto for i in 1:ndesign]

    # making a plot
    violin_plot = @pipe final_imb |> 
        stack(_, variable_name = :design) |>
        groupby(_, :design) |>
        combine(_, :value => mean => :mean) |>
        transform(_, :design => x -> categorical(x, levels = designs) => :design) |>
        scatter(_.design, zeros(ndesign), label = "",
            color = colors, 
            markersize = 2,
            ylabel = "final imbalance", 
            xrotation = 30;
            kwargs...
        )

    for i in eachindex(designs)
        x = [designs[i] for _ in 1:nsim]
        y = final_imb[:, i]
        StatsPlots.violin!(violin_plot, x, y, fillcolor = colors[i], label = "")
    end
                   
    violin_plot
end