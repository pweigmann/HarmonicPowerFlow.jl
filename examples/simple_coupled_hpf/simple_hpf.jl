# simple example of a harmonic power flow execution 
using HarmonicPowerFlow
using Plots


harmonics = [h for h in 1:2:19]
# settings
coupled_small = init_settings(true, harmonics)

# import nodes and lines from csv files
net2 = init_power_grid(
    import_nodes_from_csv("examples\\simple_coupled_hpf\\net2"), import_lines_from_csv("examples\\simple_coupled_hpf\\net2"),
    coupled_small)

# run harmonic power flow
u, err_h_final, n_iter_h = hpf(net2, coupled_small)

THD = HarmonicPowerFlow.THD(u, net2.nodes, coupled_small.harmonics)

# show barplot of harmonics at bus 4
bar(harmonics, [u[h][4, "v"] for h in harmonics], ticks = harmonics, label = "Harmonic voltages at bus 4", title = "THD = "*string(round(THD.THD_F[4]; digits= 3)))

