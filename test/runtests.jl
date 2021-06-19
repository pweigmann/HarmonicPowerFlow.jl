using Test
using HarmonicPowerFlow
using DataFrames

@testset "admittance_matrices" begin
    nodes = DataFrame(
        ID = 1:5, 
        type = ["slack", "PQ", "PQ", "PQ", "nonlinear"], 
        component = ["generator", "lin_load_1", "lin_load_2", nothing, "smps"],
        S = [1000, nothing, nothing, nothing, nothing],
        X_shunt = [0.0001, 0, 0, 0, 0],
        P = [nothing, 100, 100, 0, 250],
        Q = [nothing, 100, 100, 0, 100])
    lines = DataFrame(
        ID = 1:5,
        fromID = 1:5,
        toID = [2,3,4,5,1],
        R = [0.01, 0.02, 0.01, 0.01, 0.01],
        X = [0.01, 0.08, 0.02, 0.02, 0.02])
    settings = HarmonicPowerFlow.init_settings(true, [1, 3, 5])
    net = HarmonicPowerFlow.init_power_grid(nodes, lines, settings)
    LY = HarmonicPowerFlow.admittance_matrices(net, [1, 3, 5])
    @test size(LY[1]) == (5, 5)
    @test length(LY) == 3
end