# simple example of a harmonic power flow execution 
using HarmonicPowerFlow
using Plots


harmonics = [h for h in 1:2:21]
# settings
coupled_small = init_settings(true, harmonics)

# import nodes and lines from csv files
net2 = init_power_grid(
    import_nodes_from_csv("examples\\simple_coupled_hpf\\net2_buses.csv"), import_lines_from_csv("examples\\simple_coupled_hpf\\net2_lines.csv"),
    coupled_small)

# run harmonic power flow
u, err_h_final, n_iter_h = hpf(net2, coupled_small)

THD = HarmonicPowerFlow.THD(u)

# show barplot of harmonics at bus 4
HarmonicPowerFlow.barplot_THD(u, 4, h_max= 41)

