using Plots
using LinearAlgebra


""" Barplot of harmonic voltage magnitude at one bus """
function barplot_THD(u, nodeID; h_max=0)
    if h_max == 0
        h_max = maximum(keys(u))
    end
    h_max_idx = Int((h_max+1)/2)
    x = sort(collect(keys(u)))[1:h_max_idx]
    y = [u[h][nodeID, "v"] for h in x]

    THD = HarmonicPowerFlow.THD(u)
    bar(x, y, 
        xticks = x,
        yticks = 0:0.2:1,
        xlabel = "Harmonic Frequency", 
        ylabel = "V [pu]",
        label = "Harmonic voltages at bus 4", 
        title = "THD = "*string(round(100*THD.THD_F[nodeID]; digits= 2))*"%")
end


""" Heatmap indicating magnitude of Norton Equivalent parameters """
function heatmap_NE(nodes, settings; save=false)
    harmonics = settings.harmonics
    NE = HarmonicPowerFlow.import_Norton_Equivalents(nodes, settings)
    l = @layout [a{0.66w}  b{1h}]
    for device in collect(keys(NE))
        device = string(device)
        if settings.coupled
            p1 = heatmap(harmonics, harmonics, abs.(NE[device][2]), ticks = harmonics, title = "Y_Norton Matrix ("*device*", coupled)", yflip=true)
            p2 = heatmap(transpose(abs.(NE[device][1])), aspect_ratio=1, yflip=true, title = "I_Norton ("*device*")", yticks = nothing, border = :none)
            p = plot(p1, p2, layout = l)
        else
            p1 = heatmap(harmonics, harmonics, diagm(vec(abs.(NE[device][2]))), ticks = harmonics, title = "Y_Norton Vector ("*device*", uncoupled)", yflip=true)
            p2 = heatmap(transpose(abs.(NE[device][1])), aspect_ratio=1, yflip=true, title = "I_Norton ("*device*")", yticks = nothing, border = :none)
            p = plot(p1, p2, layout = l)
        end
        if save == "png"
            savefig(p, "heatmap_"*device)
            println("Saved heatmap to file heatmap_"*device*".png")
        elseif save == "pdf"
            savefig(p, "heatmap_"*device*".pdf")
            println("Saved heatmap to file heatmap_"*device*".pdf")
        end
        # show(p) doesn't work, for now only returned if saved to file
    end
end