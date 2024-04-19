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
    
    # setting colors
    color_scheme = ColorSchemes.tab10
    ncolors = length(color_scheme)
    colors = ncolors >= ndesign ? hcat([[color_scheme[i]] for i in 1:ndesign]...) : :auto

    # setting shapes
    shapes = [:circle :rect :star4 :diamond :star8 :utriangle :star7 :dtriangle :star5 :rtriangle];
    markers = ncolors >= ndesign ? shapes : :auto

    # making a plot 
    mm = Plots.mm
    @df op StatsPlots.plot(
        cols();
        size = (1200, 800),
        dpi = 300,
        legend = :outerright,
        legend_foreground_color = nothing,
        color = colors,
        marker = markers,
        markersize = 8,
        markerstrokewidth = 0.2, 
        xlabel = "allocation step",
        xlims = (0, nsbj+2),
        yguidefontsize = 18,
        ytickfontsize = 16,
        xguidefontsize = 18,
        xtickfontsize = 16,
        legendfontsize = 18,
        left_margin = 10mm,
        bottom_margin = 10mm,
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

    # getting number of subjects
    nsbj = size(prb, 1)

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
    mm = Plots.mm 
    hline(ρ, ls = :dash, lw = 2.0, lc = :black, label = latexify("ρ_k"))
    @df prb_df StatsPlots.plot!(
        cols();
        size = (1200, 800),
        dpi = 300,
        legend = :outerright,
        legend_foreground_color = nothing,
        xlabel = "allocation step",
        ylabel = "unconditional allocation probability",
        title = lbl, 
        xlims = (0, nsbj+2),
        ylims = (0.0, 1.0),
        xticks = :auto,
        yticks = 0:0.1:1,
        color = colors,
        marker = markers, 
        markercolor = colors,
        markersize = 8,
        markerstrokewidth = 0.5, 
        yguidefontsize = 18,
        ytickfontsize = 16,
        xguidefontsize = 18,
        xtickfontsize = 16,
        legendfontsize = 18,
        left_margin = 10mm,
        bottom_margin = 10mm,
        kwargs...
    )
end


"""Function used to visualize simulated operational characteristics.

# Call
`plot(op; kwargs...)`

# Arguments
- `op::Vector{ARP}` -- a vector of instances of `ARP`, each containing unconditional allocation probabilities.
- `kwargs` refers to the _kew words_. Here, it is possible to pass _key-value_ pairs supported by a `StatsPlots.plot` function.

# Result
- A plot of corresponding operational chacteristics vs. allocation step
"""
function plot(op::Vector{ARP}; kwargs...)
    mm = Plots.mm
    plt = [plot(op[i]) for i in eachindex(op)]
    Plots.plot(plt...; 
        size = (2400, 1800), dpi = 300,
        left_margin = 20mm, 
        bottom_margin = 20mm, 
        top_margin = 20mm, 
        right_margin = 20mm,  
        markersize = 5, 
        titlefontsize = 20,
        kwargs...)
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
    mm = Plots.mm
    @df brt_transfromed[:, Not(:design)] StatsPlots.heatmap(
        cols();
        size = (1200, 800),
        dpi = 300,
        legend = :outerright,
        xlabel = "allocation step",
        #xticks = [1; 5:5:nsbj],
        yticks = (1:ndesign, designs), 
        color = colors,
        ytickfontsize = 16,
        xguidefontsize = 18,
        xtickfontsize = 16,
        left_margin = 10mm,
        bottom_margin = 10mm,
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
        StatsPlots.scatter(_.design, zeros(ndesign), label = "",
            size = (1200, 800),
            left_margin = 10StatsPlots.mm,
            bottom_margin = 10StatsPlots.mm,
            color = colors, 
            markersize = 3,
            ylabel = "final imbalance", 
            yguidefontsize = 18,
            ytickfontsize = 16,
            xrotation = 45,
            xtickfontsize = 16;
            kwargs...
        )

    for i in eachindex(designs)
        x = [designs[i] for _ in 1:nsim]
        y = final_imb[:, i]
        StatsPlots.violin!(violin_plot, x, y, fillcolor = colors[i], label = "")
    end
                   
    violin_plot
end