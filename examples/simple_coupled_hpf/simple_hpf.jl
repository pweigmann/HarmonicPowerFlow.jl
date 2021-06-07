# simple example of a harmonic power flow execution 
using HarmonicPowerFlow

# import nodes and lines from csv files
net2 = init_power_grid(import_nodes_from_csv("examples\\simple_coupled_hpf\\net2"), import_lines_from_csv("examples\\simple_coupled_hpf\\net2"))

# settings
coupled_small = init_settings(true, [1, 3, 5])

# run harmonic power flow
u, err_h_final, n_iter_h = hpf(net2, coupled_small)