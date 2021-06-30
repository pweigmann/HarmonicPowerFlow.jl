# simple example of a harmonic power flow execution 
using HarmonicPowerFlow
using Plots
using LinearAlgebra


harmonics = [h for h in 1:2:5]
# settings
coupled = true
settings1 = init_settings(coupled, harmonics)

# import nodes and lines from csv files
net3 = init_power_grid(
    import_nodes_from_csv("examples\\multiple_pv\\net3"), import_lines_from_csv("examples\\multiple_pv\\net3"),
    settings1)

# run harmonic power flow
u, err_h_final, n_iter_h = hpf(net2, settings1)

THD = HarmonicPowerFlow.THD(u, net2.nodes, settings1.harmonics)

# show barplot of harmonics at bus 4
bar(harmonics, [u[h][4, "v"] for h in harmonics], ticks = harmonics, label = "Harmonic voltages at bus 4", title = "THD = "*string(round(THD.THD_F[4]; digits= 3)))

# show barplot of harmonics at bus 3
bar(harmonics, [u[h][3, "v"] for h in harmonics], ticks = harmonics, label = "Harmonic voltages at bus 3", title = "THD = "*string(round(THD.THD_F[3]; digits= 3)))

# Analyze Norton parameters
NE = HarmonicPowerFlow.import_Norton_Equivalents(net2.nodes, settings1)
l = @layout [a{0.66w}  b{1h}]
if coupled
    p1 = heatmap(harmonics, harmonics, abs.(NE["SMPS"][2]), ticks = harmonics, title = "Y_Norton Matrix (coupled)", yflip=true)
    p2 = heatmap(transpose(abs.(NE["SMPS"][1])), aspect_ratio=1, yflip=true, title = "I_Norton", yticks = nothing, border = :none)
    plot(p1, p2, layout = l)
else
    Y_N = heatmap(harmonics, harmonics, diagm(vec(abs.(NE["SMPS"][2]))), ticks = harmonics, title = "Y_Norton Vector (uncoupled)", yflip=true)
end

[[h, u[h]] for h in harmonics]