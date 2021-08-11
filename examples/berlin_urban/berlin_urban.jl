# simple example of a harmonic power flow execution 
using HarmonicPowerFlow
using LinearAlgebra
using StatsPlots
using DataFrames
using CSV

harmonics = [h for h in 1:2:101]

# settings
settings_c = init_settings(true, harmonics, base_voltage=400)
settings_uc = init_settings(false, harmonics, base_voltage=400)

# import nodes and lines from csv files
urban_12 = init_power_grid(
    import_nodes_from_csv("examples\\berlin_urban\\berlin_urban_buses_12.csv"), 
    import_lines_from_csv("examples\\berlin_urban\\berlin_urban_lines_12.csv"),
    settings_c)
urban_36 = init_power_grid(
    import_nodes_from_csv("examples\\berlin_urban\\berlin_urban_buses_36.csv"), 
    import_lines_from_csv("examples\\berlin_urban\\berlin_urban_lines_36.csv"),
    settings_c)
urban_72 = init_power_grid(
    import_nodes_from_csv("examples\\berlin_urban\\berlin_urban_buses_72.csv"), 
    import_lines_from_csv("examples\\berlin_urban\\berlin_urban_lines_72.csv"),
    settings_c)
urban_mix = init_power_grid(
    import_nodes_from_csv("examples\\berlin_urban\\berlin_urban_buses_72_mix.csv"), 
    import_lines_from_csv("examples\\berlin_urban\\berlin_urban_lines_72.csv"),
    settings_c)

u_c_12,  err_h_c_12, n_iter_h = hpf(urban_12, settings_c)
u_uc_12,  err_h_c_12, n_iter_h = hpf(urban_12, settings_uc)
u_c_36,  err_h_c_36, n_iter_h = hpf(urban_36, settings_c)
u_uc_36, err_h_uc_36, n_iter_h = hpf(urban_36, settings_uc)
u_c_72,  err_h_c_72, n_iter_h = hpf(urban_72, settings_c)
u_uc_72, err_h_uc_72, n_iter_h = hpf(urban_72, settings_uc)
u_c_mix,  err_h_c_mix, n_iter_h = hpf(urban_mix, settings_c)
u_uc_mix, err_h_uc_mix, n_iter_h = hpf(urban_mix, settings_uc)