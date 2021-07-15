using DataFrames

"""
    THD(u)

Calculate the Total Harmonic Distortion at all nodes.

THD is a measure for the amount of voltage distortion present at a node. 
This function calculates two alternative definitions 
    - THD_F: relative to fundamental, more commonly used and 
    - THD_R: RMS, relative to all frequencies.
"""
function THD(u)
    harmonics = sort(collect(keys(u)))
    THD = DataFrame(THD_F = zeros(size(u[1], 1)), 
                    THD_R = zeros(size(u[1], 1)))
    for ID in 1:size(u[1], 1)
        THD[ID, "THD_F"] = sqrt(sum([u[h].v[ID]^2 for h in harmonics[2:end]]))./u[1].v[ID]
        THD[ID, "THD_R"] = sqrt(sum([u[h].v[ID]^2 for h in harmonics[2:end]]))./sqrt(sum([u[h].v[ID]^2 for h in harmonics]))
    end
    return THD
end


"""
    EN50160(u, nodeID)

Check if harmonic limits according to EN 50160 and EN 61000 are met at node."""
function limits(u, nodeID)
    h = [3, 5, 7, 9, 11, 13, 15]
    v = [u[i].v[nodeID] for i in h]
    limit_EN50160 = [0.05, 0.06, 0.05, 0.015, 0.035, 0.03, 0.005]
    limit_EN61000 = [0.05, 0.06, 0.05, 0.015, 0.035, 0.03, 0.005]
    met_EN50160 = collect(v .< limit_EN50160)
    met_EN61000 = collect(v .< limit_EN61000)
    
    if all(met_EN50160) && all(met_EN61000)
        println("Harmonic limits met at node "*string(nodeID))
    else
        println("Harmonic limits not met at node "*string(nodeID)*"! Check DataFrame for details.")
    end

    THD = HarmonicPowerFlow.THD(u).THD_F[nodeID]

    if THD < 0.08
        println("THD at node "*string(nodeID)*" within limits (< 8%).")
    else
        println("THD at node "*string(nodeID)*" exceeds limits (> 8%)!")
    end

    DataFrame(v=v, within_EN50160 = met_EN50160, within_EN61000=met_EN61000)
end


# function TDD()

# function PF()