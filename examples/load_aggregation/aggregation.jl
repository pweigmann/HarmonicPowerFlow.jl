# simple example of a harmonic power flow execution 
using HarmonicPowerFlow
using StatsPlots


harmonics = [h for h in 1:2:101]

# settings
settings_c = init_settings(true, harmonics)
settings_uc= init_settings(false, harmonics)

# import nodes and lines from csv files
net_agg1 = init_power_grid(
    import_nodes_from_csv("examples\\load_aggregation\\aggregated_buses_case1.csv"), import_lines_from_csv("examples\\load_aggregation\\aggregated_lines_case1.csv"),
    settings_c)

net_agg2 = init_power_grid(
    import_nodes_from_csv("examples\\load_aggregation\\aggregated_buses_case2.csv"), import_lines_from_csv("examples\\load_aggregation\\aggregated_lines_case2.csv"),
    settings_c)

# run harmonic power flow
u1_c, err_h_final1, n_iter_h1 = hpf(net_agg1, settings_c)
u2_c, err_h_final2, n_iter_h2 = hpf(net_agg2, settings_c)
u1_uc, err_h_final1, n_iter_h1 = hpf(net_agg1, settings_uc)
u2_uc, err_h_final2, n_iter_h2 = hpf(net_agg2, settings_uc)

THD1_c = HarmonicPowerFlow.THD(u1_c)
THD2_c = HarmonicPowerFlow.THD(u2_c)
THD1_uc = HarmonicPowerFlow.THD(u1_uc)
THD2_uc = HarmonicPowerFlow.THD(u2_uc)

# show harmonics at bus 2 (worst case)
HarmonicPowerFlow.barplot_THD(u1_c, 1, h_max= 41)
HarmonicPowerFlow.barplot_THD(u2_c, 1, h_max= 41)
HarmonicPowerFlow.barplot_THD(u1_uc, 1, h_max= 41)
HarmonicPowerFlow.barplot_THD(u2_uc, 1, h_max= 41)

# filename_results = "agg_test_sc2_bus2"
# Plots.savefig(p, filename_results)