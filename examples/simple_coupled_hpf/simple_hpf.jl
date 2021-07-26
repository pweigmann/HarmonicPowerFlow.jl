# simple example of a harmonic power flow execution 
using HarmonicPowerFlow
using StatsPlots


harmonics = [h for h in 1:2:21]
# settings
settings_c = init_settings(true, harmonics)
settings_uc = init_settings(false, harmonics)

# import nodes and lines from csv files
grid = init_power_grid(
    import_nodes_from_csv("examples\\simple_coupled_hpf\\net2_buses.csv"), 
    import_lines_from_csv("examples\\simple_coupled_hpf\\net2_lines.csv"),
    settings_uc)

# run harmonic power flow
u, err_h_final, n_iter_h = hpf(grid, settings_uc)

THD = HarmonicPowerFlow.THD(u)

# show barplot of harmonics at bus 4
HarmonicPowerFlow.barplot_THD(u, 4, h_max= 21)

