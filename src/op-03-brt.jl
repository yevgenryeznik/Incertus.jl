"""Function calculates balance-randomness trade-off vs. allocation step.

# Call
`calc_brt(sr)`

# Arguments
- `sr::SimulatedRandomization`: an instance of `SimulatedRandomization`, an object, 
representing simulation output.

# Result
- A `Vector` of _balance-randomness trade-off measurement_ values summarized via simulations.
"""
function calc_brt(sr::SimulatedRandomization)
    # extracting information about target allocation
    target = sr.target

    # extracting informaiton about simulation parameters
    nsbj, ntrt, nsim = size(sr.trt)

    # measure of imbalance vs. allocation step
    cummean_loss = calc_cummean_loss(sr)

    # measure of randomness vs. allocation step
    fi = calc_fi(sr)
    
    if (ntrt == 2) && (target[1] == target[2])
        # balance-randomness trade-off vs. allocation step
        G = [sqrt(item1^2 + item2^2) for (item1, item2) in zip(cummean_loss, fi)]
    else
        cummean_loss_crd = Float64[]
        fi_crd = Float64[]

        cummean_loss_maxent = Float64[]
        fi_maxent = Float64[]

        if (sr.label != "CRD")
           println("simulating CRD needed for evaluating balance-randomness trade-off...")
           crd_simulated = simulate(CRD(target), nsbj, nsim) 

           cummean_loss_crd = calc_cummean_loss(crd_simulated)
           fi_crd = calc_fi(crd_simulated)
        else
            cummean_loss_crd = cummean_loss
            fi_crd = fi 
        end    
        if (sr.label != "MaxEnt(1)")
            println("simulating MaxEnt(1) needed for evaluating balance-randomness trade-off...")
            maxent_simulated = simulate(MaxEnt(target, 1), nsbj, nsim)

            cummean_loss_maxent = calc_cummean_loss(maxent_simulated)
            fi_maxent = calc_fi(maxent_simulated) 
        else
            cummean_loss_maxent = cummean_loss
            fi_maxent = fi 
        end    
        UI = (cummean_loss - cummean_loss_maxent) ./ (cummean_loss_crd - cummean_loss_maxent)
        UR = (fi - fi_crd) ./ (fi_maxent - fi_crd)

        G = [sqrt((UI[i]^2 + UR[i]^2)) for i in 1:nsbj]
    end

    return G
end


"""Function calculates balance-randomness trade-off vs. allocation step.

# Call
`calc_brt(sr)`

# Arguments
- `sr::Vector{SimulatedRandomization}`: a vector of instances of `SimulatedRandomization`, representing simulation output.

# Result
- A `DataFrame` of _balance-randomness trade-off measurements_' values summarized via simulations.
"""
function calc_brt(sr::Vector{SimulatedRandomization})
    # extracting information about target allocation
    target = sr[1].target

    # extracting informaiton about simulation parameters
    nsbj, ntrt, nsim = size(sr[1].trt)
            
    if (ntrt == 2) && (target[1] == target[2])
        brt = hcat([calc_brt(item) for item in sr]...)
        rnd_labels = [item.label for item in sr]

        return DataFrame(brt, rnd_labels)
    else
        rnd_labels = [item.label for item in sr]

        # measure of imbalance vs. allocation step
        cummean_loss = calc_cummean_loss(sr)
    
        # measure of randomness vs. allocation step
        fi = calc_fi(sr)

        cummean_loss_crd = Float64[]
        fi_crd = Float64[]

        cummean_loss_maxent = Float64[]
        fi_maxent = Float64[]
        
        if ("CRD" in rnd_labels)
            cummean_loss_crd = cummean_loss[:, "CRD"]
            fi_crd = fi[:, "CRD"]
        else
            println("simulating CRD needed for evaluating balance-randomness trade-off...")
            crd_simulated = simulate(CRD(target), nsbj, nsim)

            cummean_loss_crd = calc_cummean_loss(crd_simulated)
            fi_crd = calc_fi(crd_simulated)
        end

        if ("MaxEnt(1)" in rnd_labels)
            cummean_loss_maxent = cummean_loss[:, "MaxEnt(1)"]
            fi_maxent = fi[:, "MaxEnt(1)"]
        else
            println("simulating MaxEnt(1) needed for evaluating balance-randomness trade-off...")
            maxent_simulated = simulate(MaxEnt(target, 1), nsbj, nsim)
             
            cummean_loss_maxent = calc_cummean_loss(maxent_simulated)
            fi_maxent = calc_fi(maxent_simulated)
        end
        UI = hcat([(col - cummean_loss_maxent) ./ (cummean_loss_crd - cummean_loss_maxent) for col in eachcol(cummean_loss)]...)
        UR = hcat([(col - fi_crd) ./ (fi_maxent - fi_crd) for col in eachcol(fi)]...)

        G = sqrt.((UI.^2 + UR.^2))

        return DataFrame(G, rnd_labels)
    end
end