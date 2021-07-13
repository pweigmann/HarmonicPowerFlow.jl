using DataFrames
using CSV
using SparseArrays

"""
    admittance_matrices(nodes, lines, harmonics)

Build the nodal admittance matrices (admittance laplacian) for all harmonics. Admittance scales with frequency (X_h = X_f * h).

Returns: dictionary of DataFrames LY[h]
"""
function admittance_matrices(net, harmonics)
    LY = Dict()
    for h in harmonics
        LY[h] = spzeros(ComplexF64, net.n, net.n)
        # non-diagonal elements
        for line in eachrow(net.lines)
            LY[h][line.fromID, line.toID] = -1/(line.R + 1im*line.X*h)
            # nodal admittance matrix is assumed to be symmetric
            LY[h][line.toID, line.fromID] = -1/(line.R + 1im*line.X*h)
        end
        # diagonal elements
        for i in 1:net.n
            # bus shunt admittances only for harmonic frequencies
            if net.nodes.X_sh[i] > 0 && h != 1
                LY[h][i, i] = -sum(LY[h][i, :]) + 1/(1im*net.nodes.X_sh[i]*h)
            else
                LY[h][i, i] = -sum(LY[h][i, :])
            end
            # Adding shunt admittances for each pi-model line connected to bus i
            for line in eachrow(net.lines)
                # check if line is connected to bus i
                if line.fromID == i || line.toID == i
                    LY[h][i, i] = LY[h][i, i] + 
                    (line.G + 1im*h*line.B)/2
                end
            end
        end
    end
    return LY
end


function init_voltages(nodes, settings)
    u = Dict()
    for h in settings.harmonics
        if h == 1
            u[h] = DataFrame(
                v = ones(size(nodes, 1))*settings.v_f,
                ϕ = ones(size(nodes, 1))*settings.ϕ_f
            )  
        else
            u[h] = DataFrame(
                v = ones(size(nodes, 1))*settings.v_h,
                ϕ = ones(size(nodes, 1))*settings.ϕ_h
            )
        end
    end
    return u
end


"""
    fund_state_vec(u)

Take the voltage DataFrame and return a complex vector.

Note: sorting of values. For the fundamental state vector voltage phase comes first, then magnitude.
"""
function fund_state_vec(net, u)
    xϕ = u[1].ϕ[2:end]
    xv = u[1].v[(net.c+1):end]
    vcat(xϕ, xv) 
end


function fund_mismatch(net, u, LY)
    LY_1 = LY[1]
    u_1 = u[1].v .* exp.(1im*u[1].ϕ)
    s = (net.nodes.P + 1im*net.nodes.Q)
    mismatch = u_1 .* conj(LY_1*u_1) + s
    f = vcat(real(mismatch[2:end]), imag(mismatch[(net.c+1):end]))
    err = maximum(abs.(f))
    return f, err
end


function fund_jacobian(net, u, LY) 
    u_1 = u[1].v .* exp.(1im*u[1].ϕ)
    i_diag = spdiagm(sparse(LY[1] * u_1))
    u_1_diag = spdiagm(u_1)
    u_1_diag_norm = spdiagm(u_1./abs.(u_1))

    dSdϕ = 1im*u_1_diag*conj(i_diag - LY[1]*u_1_diag)
    dSdv = u_1_diag_norm*conj(i_diag) + u_1_diag*conj(LY[1]*u_1_diag_norm)

    # divide sub-matrices into real and imag part, cut off slack
    dPdϕ = real(dSdϕ[2:end, 2:end])
    dPdv = real(dSdv[2:end, (net.c+1):end])
    dQdϕ = imag(dSdϕ[(net.c+1):end, 2:end])
    dQdv = imag(dSdv[(net.c+1):end, (net.c+1):end])

    vcat(hcat(dPdϕ, dPdv), 
         hcat(dQdϕ, dQdv))
end


function update_state_vec!(J, x, f)
    x - J\f  # Newton-Raphson iteration
    # -> could try using NLsolve.jl
end


function update_fund_voltages!(net, u, x)
    u[1].ϕ[2:end] = x[1:(net.n-1)]
    u[1].v[(net.c+1):end] = x[net.n:end]
    return u
end


function pf(net, settings, LY, plt_convergence = false)
    u = init_voltages(net.nodes, settings)
    x_f = fund_state_vec(net, u)
    f_f, err_f = fund_mismatch(net, u, LY)

    n_iter_f = 0
    err_f_t = Dict()
    while err_f > settings.thresh_f && n_iter_f <= settings.max_iter_f
        J_f = fund_jacobian(net, u, LY)
        x_f = update_state_vec!(J_f, x_f, f_f)
        u = update_fund_voltages!(net, u, x_f)
        f_f, err_f = fund_mismatch(net, u, LY)
        err_f_t[n_iter_f] = err_f
        n_iter_f += 1
    end

    println(u[1])
    if n_iter_f < settings.max_iter_f
        println("Fundamental power flow converged after ", n_iter_f, 
              " iterations (err < ", settings.thresh_f, ").")
    elseif n_iter_f == settings.max_iter_f
        println("Maximum of ", n_iter_f, " iterations reached.")
    end
    return u
end


# Harmonic Power Flow functions
"""Import Norton Equivalent parameters for all nonlinear devices in "nodes" from folder"""
function import_Norton_Equivalents(nodes, settings, folder_path="devices\\")
    NE = Dict()
    nl_components = unique(nodes[nodes.type .== "nonlinear", "component"])
    for device in nl_components
        NE_df = CSV.read(folder_path * device * "_NE.csv", DataFrame)
        # transform to Complex type, enough to strip first paranthesis for successful parse
        vals = mapcols!(col -> parse.(ComplexF64, strip.(col, ['('])), NE_df[:, 3:end])
        NE_device = hcat(NE_df[:,1:2], vals)
        # filter columns for considered harmonics
        NE_device = NE_device[:, Between(begin, string(maximum(settings.harmonics)*settings.base_frequency))]
        # --> CHECK if this works as intended

        # change to pu system and choose if coupled 
        if settings.coupled
            I_N = Array(NE_device[NE_device.Parameter .== "I_N_c", 3:end])/settings.base_current
            LY_N_full = NE_device[NE_device.Parameter .== "Y_N_c", 2:end]
            LY_N = Array(LY_N_full[LY_N_full.Frequency .<= maximum(settings.harmonics)*settings.base_frequency, 2:end])/settings.base_admittance
        else
            I_N = Array(NE_device[NE_device.Parameter .== "I_N_uc", 3:end])/settings.base_current
            LY_N = Array(NE_device[NE_device.Parameter .== "Y_N_uc", 3:end])/settings.base_admittance
        end    
        NE[device] = [I_N, LY_N]
    end
    return NE
end


"""calculate the harmonic current injections at one node"""
function current_injections(nodes, nodeID, u, NE, harmonics)
    component = nodes[nodes.ID .== nodeID, "component"][1]
    I_N, LY_N = NE[component]
    # u as dict of dfs makes building this vector a bit complicated
    u_h = vcat([u[h][nodeID, "v"] .* exp.(1im*u[h][nodeID, "ϕ"]) for h in harmonics]...)
    # coupled: Y_N is a matrix, uncoupled: vector
    if size(LY_N)[1] > 1  # coupled case
        i_inj = vec(I_N) - vec(LY_N*u_h)
    else  # uncoupled case
        i_inj = vec(I_N) - spdiagm(vec(LY_N))*u_h
    end
    return i_inj
end

function current_balance(net, settings, u, LY, NE)
    # fundamental admittance matrix for nonlinear nodes
    LY_1_nl = LY[1][net.m:end,:]
    u_1 = u[1].v .* exp.(1im*u[1].ϕ)
    # fundamental line currents at nonlinear nodes
    dI_1 = LY_1_nl * u_1
    # harmonic admittance matrices as diagonal block matrix
    LY_h = blockdiag([LY[h] for h in settings.harmonics[2:end]]...)
    u_h = vcat([u[h][:, "v"] .* exp.(1im*u[h][:, "ϕ"]) for h in settings.harmonics[2:end]]...)
    dI_h = LY_h * u_h 

    # subtract the injected currents at each nonlinear node i
    for i in net.m:net.n
        i_inj = current_injections(net.nodes, net.nodes.ID[i], u, NE, settings.harmonics)
        dI_1[i-net.m+1] += i_inj[1]  # add injections at fundamental frequency...
        # ... and at all harmonic frequencies
        for p in 0:(settings.K-1)
            dI_h[p*net.n + i] += i_inj[p+2]
        end
    end
    vcat(dI_1, dI_h)
end


function harmonic_mismatch(net, settings, u, LY, NE)
    # fundamental power mismatch at linear buses except slack
    s = net.nodes.P[2:(net.m-1)] + 1im*net.nodes.Q[2:(net.m-1)]
    u_i = u[1][2:(net.m-1), "v"].*exp.(1im*u[1][2:(net.m-1), "ϕ"])
    u_j = u[1][:, "v"].*exp.(1im*u[1][:, "ϕ"])
    LY_ij = LY[1][2:(net.m-1), :]
    # power balance
    sl = u_i.*conj(LY_ij*u_j)  
    ds = s + sl  
    di = current_balance(net, settings, u, LY, NE) 
    # harmonic mismatch vector
    f_c = vcat(ds, di) 
    f = vcat(real(f_c), imag(f_c[(net.c):end]))  # crop PV buses
    err_h = maximum(abs.(f))  
    return f, err_h
end


function harmonic_state_vec(net, u, harmonics)
    # fundamental voltages, crop slack and magnitude for PV buses
    xϕ = u[1].ϕ[2:end]
    xv = u[1].v[(net.c+1):end]
    # all voltages and angles at all harmonic frequencies 
    for h in harmonics[2:end]
        xϕ = vcat(xϕ, u[h].ϕ)
        xv = vcat(xv, u[h].v)
    end
    vcat(xϕ, xv)  # note: phase first
end   


function build_harmonic_jacobian(net, settings, u, LY, NE)
    m = net.m
    n = net.n
    c = net.c
    K = settings.K

    u_vec = vcat([u[h].v .* exp.(1im*u[h].ϕ) for h in settings.harmonics]...)
    v_vec = vcat([u[h].v for h in settings.harmonics]...)
    u_diag = spdiagm(u_vec)
    u_norm = u_vec./v_vec
    u_norm_diag = spdiagm(u_norm)
    LY_diag = blockdiag([LY[h] for h in settings.harmonics]...)

    # construct Jacobian sub-matrices
    dIdv = LY_diag*u_norm_diag
    dIdϕ = 1im*LY_diag*u_diag

    # indices of first nonlinear bus at each harmonic
    nl_idx_start = m:n:n*(K+1)
    nl_idx_all = vcat([nl:(nl+n-m) for nl in nl_idx_start]...)
    u_nl = u_vec[nl_idx_all]
    v_nl = v_vec[nl_idx_all]
    u_nl_norm = u_nl./v_nl

    if settings.coupled
        # iterating through blocks vertically
        for h in 0:K
            # ... and horizontally
            for p in 0:K
                # iterating through nonlinear buses
                for i in m:n
                    # within NE "[2]" points to LY_N
                    LY_N = NE[net.nodes.component[i]][2]
                    # subtract derived current injections at respective idx
                    dIdv[h*n+i, p*n+i] -= LY_N[h+1, p+1]*u_nl_norm[(i-m+1)+p*(n-m+1)]
                    dIdϕ[h*n+i, p*n+i] -= 1im*LY_N[h+1, p+1]*u_nl[(i-m+1)+p*(n-m+1)]
                end
            end
        end
    else
        # iterating through blocks diagonally (p=h)
        for h in 0:K
            for i in m:n
                # LY_N is one-dimensional in the uncoupled case
                LY_N = NE[net.nodes.component[i]][2]
                dIdv[h*n+i, h*n+i] -= LY_N[h+1]*u_nl_norm[(i-m+1)+h*(n-m+1)]
                dIdϕ[h*n+i, h*n+i] -= 1im*LY_N[h+1]*u_nl[(i-m+1)+h*(n-m+1)]
            end
        end
    end

    dIdϕ = dIdϕ[m:end, 2:end]
    dIdv = dIdv[m:end, (c+1):end]  
      
    LY_1 = LY[1]
    u_1 = u[1].v .* exp.(1im*u[1].ϕ)
    i_diag = spdiagm(LY_1*u_1)
    u_diag = spdiagm(u_1)
    u_diag_norm = spdiagm(u_1./abs.(u_1))

    dS1dϕ1 = 1im*u_diag*(conj(i_diag - LY_1*u_diag))
    dS1dv1 = u_diag_norm*conj(i_diag) + u_diag*conj(LY_1*u_diag_norm)
    
    dSdϕ = hcat(dS1dϕ1, zeros(n, n*K))
    dSdv = hcat(dS1dv1, zeros(n, n*K))
    
    dPdϕ = real(dSdϕ[2:(m-1), 2:end])
    dPdv = real(dSdv[2:(m-1), (c+1):end])
    dQdϕ = imag(dSdϕ[(c+1):(m-1), 2:end])
    dQdv = imag(dSdv[(c+1):(m-1), (c+1):end])

    # combine all sub-matrices and return complete Jacobian
    vcat(hcat(dPdϕ, dPdv),
         hcat(real(dIdϕ), real(dIdv)),
         hcat(dQdϕ, dQdv),
         hcat(imag(dIdϕ), imag(dIdv)))
end


function update_harmonic_voltages!(net, u, x, settings)
    n = net.n
    c = net.c

    # slice x (almost) in half to separate voltage magnitude and phase
    xϕ = x[1:(n*settings.K1-1)]
    xv = x[n*settings.K1:end]

    # xϕ = xϕ .% (2*π)  # ensure phase smaller 2π, not really necessary

    # update fundamental voltages
    u[1].ϕ[2:end] = xϕ[1:(n-1)]
    u[1].v[(c+1):end] = xv[1:(n-c)]

    xϕ_h = xϕ[n:end]
    xv_h = xv[(n-c+1):end]
    
    for h in settings.harmonics[2:end]
        # enumerate harmonics for easier indexing
        i = findall(settings.harmonics[2:end] .== h)[1] - 1  
        u[h].ϕ = xϕ_h[(i*n+1):(i*n+n)]
        u[h].v = xv_h[(i*n+1):(i*n+n)]
    end
    return u
end


"""
    hpf(net, settings)

Solve the harmonic power flow problem for a given power grid and settings."""
function hpf(net, settings)
    LY = admittance_matrices(net, settings.harmonics)
    u = pf(net, settings, LY)
    NE = import_Norton_Equivalents(net.nodes, settings)
    f, err_h = harmonic_mismatch(net, settings, u, LY, NE)
    x = harmonic_state_vec(net, u, settings.harmonics)
    n_iter_h = 0
    err_h_t = Dict()
    while err_h > settings.thresh_h && n_iter_h < settings.max_iter_h
        J = build_harmonic_jacobian(net, settings, u, LY, NE)
        x = update_state_vec!(J, x, f)
        u = update_harmonic_voltages!(net, u, x, settings)  
        f, err_h = harmonic_mismatch(net, settings, u, LY, NE)
        err_h_t[n_iter_h] = err_h
        n_iter_h += 1
    end

    # getting rid of negative voltage magnitudes:
    for h in settings.harmonics
        u[h].ϕ[u[h].v .< 0] .+= π
        u[h].ϕ .= u[h].ϕ .% (2*π)
        u[h].v[u[h].v .< 0] = -u[h].v[u[h].v .< 0]
    end

    if n_iter_h < settings.max_iter_h
        println("Harmonic power flow converged after ", n_iter_h,
                " iterations (err < ", settings.thresh_h, ").")
    elseif n_iter_h == settings.max_iter_h
        println("Maximum of ", n_iter_h, " iterations reached. Harmonic power flow did not converge.")
    end
    return u, err_h, n_iter_h
end