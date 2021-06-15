function get_THD(u, nodes, harmonics)
    THD = DataFrame(THD_F = zeros(size(nodes)[1]), 
                    THD_R = zeros(size(nodes)[1]))
    for ID in nodes.ID
        THD[ID, "THD_F"] = sqrt(sum([u[h].v[ID]^2 for h in harmonics[2:end]]))./u[1].v[ID]
        THD[ID, "THD_R"] = sqrt(sum([u[h].v[ID]^2 for h in harmonics[2:end]]))./sqrt(sum([u[h].v[ID]^2 for h in harmonics]))
    end
end