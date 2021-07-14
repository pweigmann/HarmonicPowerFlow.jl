using Plots

""" Barplot of harmonic voltage magnitude at one bus """
function barplot_THD(u, node; h_max=0)
    if h_max == 0
        h_max = maximum(keys(u))
    end
    h_max_idx = Int((h_max+1)/2)
    x = sort(collect(keys(u)))[1:h_max_idx]
    y = [u[h][node, "v"] for h in x]

    THD = HarmonicPowerFlow.THD(u)
    bar(x, y, 
        ticks = x, 
        label = "Harmonic voltages at bus 4", 
        title = "THD = "*string(round(100*THD.THD_F[node]; digits= 2))*"%")
end