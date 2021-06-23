using DataFrames
using CSV

# Structs
struct PowerGrid
    nodes::DataFrame
    lines::DataFrame
    m::Int
    n::Int
end


"""
    init_power_grid(nodes, lines)

Create a PowerGrid struct from nodes and lines DataFrames and converts all quantities to the p.u. system. 
"""
function init_power_grid(nodes::DataFrame, lines::DataFrame, settings)
    m = minimum(nodes[nodes[:,"type"] .== "nonlinear",:].ID)  # number of linear nodes
    n = size(nodes, 1)  # total number of nodes

    # convert quantities to p.u. system
    nodes.S = nodes.S./settings.base_power
    nodes.P = nodes.P./settings.base_power
    nodes.Q = nodes.Q./settings.base_power
    nodes.X_sh = nodes.X_sh./settings.base_impedance

    lines.R = lines.R./settings.base_impedance
    lines.X = lines.X./settings.base_impedance
    lines.G = lines.G./settings.base_admittance
    lines.B = lines.B./settings.base_admittance

    return PowerGrid(nodes, lines, m, n)
end


function import_nodes_from_csv(filename)
    df = CSV.read(filename * "_buses.csv", DataFrame)
end


function import_lines_from_csv(filename)
    CSV.read(filename * "_lines.csv", DataFrame)
end

"""
   create_nodes_manually()

Manually create a nodes DataFrame. 

Note that linear nodes are added first, then nonlinear ones. First node contains the slack bus, which is also the only node that currently is allowed to have a shunt admittance."""
function create_nodes_manually()
    DataFrame(
        ID = 1:5, 
        type = ["slack", "PQ", "PQ", "PQ", "nonlinear"], 
        component = ["generator", "lin_load_1", "lin_load_2", nothing, "smps"],
        S = [1000, nothing, nothing, nothing, nothing],
        P = [nothing, 100, 100, 0, 250],
        Q = [nothing, 100, 100, 0, 100],
        X_sh = [0.0001, 0, 0, 0, 0])
end


function create_lines_manually()
    DataFrame(
        ID = 1:5,
        fromID = 1:5,
        toID = [2,3,4,5,1],
        R = [0.01, 0.02, 0.01, 0.01, 0.01],
        X = [0.01, 0.08, 0.02, 0.02, 0.02],
        G = [0, 0, 0, 0, 0],
        B = [0, 0, 0, 0, 0])
end