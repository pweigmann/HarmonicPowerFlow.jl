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


# function TDD()

# function PF()