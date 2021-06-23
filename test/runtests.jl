using Test
using HarmonicPowerFlow
using DataFrames

@testset "admittance_matrices" begin
    nodes = DataFrame(
        ID = 1:5, 
        type = ["slack", "PQ", "PQ", "PQ", "nonlinear"], 
        component = ["generator", "lin_load_1", "lin_load_2", "lin_load_3", "smps"],
        S = [1000, 0, 0, 0, 0],
        P = [0, 100, 100, 0, 250],
        Q = [0, 100, 100, 0, 100],
        X_sh = [0.005, 0, 0, 0, 0])
    lines = DataFrame(
        ID = 1:5,
        fromID = 1:5,
        toID = [2,3,4,5,1],
        R = [0.5, 1, 0.5, 0.5, 0.5],
        X = [0.5, 4, 1, 1, 1])
    settings = HarmonicPowerFlow.init_settings(true, [1, 3, 5])
    net = HarmonicPowerFlow.init_power_grid(nodes, lines, settings)
    LY = HarmonicPowerFlow.admittance_matrices(net, [1, 3, 5])
    @test size(LY[1]) == (5, 5)
    @test length(LY) == 3
end