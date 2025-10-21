using Revise
using CSV, DataFrames
using Plots

includet("../../src/axion_red_vec/pnjl_vec.jl")


function data_Gv()

    T = 10.0

    #rhos = 0.1:0.1:13.0
    rhos = 3.00:-0.01:0.01
    theta = 0.0
    X0 = [1.8, 1.8, 2.2, -0.001, -0.001, -0.001, 0.00383, 0.00383, 1156/hc, -9.38/hc]

    #X0 = [0.43385,0.52539,2.03942,0.95168,1.02742,1.16899,0.00108,0.00123,766.34909/hc,-24.26929/hc] #theta = pi
   
    #X0 = [1.70734,1.69475,2.25882,-0.11777,-0.07,-0.07907,0.00143,0.0017,1034.32725/hc,-8.12306/hc] # theta = pi/3
    Rv = 0.0
    Gv = Rv * G_f

    unit_ratio = 2.6115e-4 # fm^-4 to km^-2
    ints = get_nodes(800)
    ints2 = ints[2]
    data = zeros(length(rhos), 11)
    EOS = zeros(length(rhos), 10)
    data[:, 1] = rhos
    EOS[:, 1] = rhos
    #P0 = ((0.42544232e3)/hc)^4
    P0 = 21.0702317062811
    for (i, rho) in enumerate(rhos)
    
        X0 = Trho(X0, T/hc ,theta, rho, Gv, ints)
        data[i, 2:end] = X0
        data[i, end-1:end] = X0[end-1:end] * hc
        phi = X0[1:6]
        Phi1 = X0[7]
        Phi2 = X0[8]
        mu_B = X0[9] * hc
        mu_Q = X0[10] * hc
        P = - Omega(phi, Phi1, Phi2, mu_B/hc, mu_Q/hc, T/hc, theta, rho, Gv, ints) - P0
        rho_eff = rho_mu_eff(mu_B/hc, mu_Q/hc, T/hc, theta, X0[1:8], rho, Gv, ints2)

        S = - dOmgea_dT(phi * 1.0001, Phi1 * 1.0001, Phi2 * 1.0001, mu_B/hc, mu_Q/hc, T/hc, theta, rho, Gv, ints)
        E = - P + T/hc * S + rho_eff[1]
        EOS[i, 2] = P * hc
        EOS[i, 3] = E * hc
        EOS[i, 4] = S * unit_ratio
        EOS[i, 5:end] .= rho_eff[1:end]
    end

    df = DataFrame(data, [:rho, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2, :mu_B, :mu_Q])
    df2 = DataFrame(EOS, [:rho, :P, :E, :S, :rho_eff, :rho_u, :rho_d, :rho_s, :rho_e, :rho_mu])
    
    outpath1 = "../../data/axion_red_vec/Rv=$(Rv)_T=$T.dat"
    outpath2 = "../../data/axion_red_vec/EOS_Rv=$(Rv)_T=$T.dat"
    CSV.write("Rv=$(Rv)_T=$T.dat", df)
    CSV.write("EOS_Rv=$(Rv)_T=$T.dat", df2)
end



function find_P0()
    T = 5.0
    mu_B = 0.0
    dθ   = 0.01
    ϵ    = 1e-6  # 避免端点重合
    θ1   = range(0.0,     step=dθ, stop=pi - ϵ)
    θ2   = range(pi,      step=dθ, stop=3pi - ϵ)
    θ3   = range(3pi,     step=dθ, stop=4pi)
    lens = length(θ1) + length(θ2) + length(θ3)


    mu_Q = -10.0


    rho = 0.1
    X0 = [1.8, 1.8, 2.2, -0.001, -0.001, -0.001, 0.00383, 0.00383, rho, mu_Q/hc]
    Rv = 0.600
    Gv = Rv * G_f
 
    X_get = zeros(lens, length(X0) + 2)
  
    X_get[:, 1] .= vcat(collect(θ1), collect(θ2), collect(θ3))

    ints = get_nodes(256)

     # 段 1: θ ∈ [0, π)
    println("Calculating segment 1...")
    for (i, θ) in enumerate(θ1)
        X0 = Tmu(X0, mu_B/hc, T/hc, θ, Gv, ints)
        orders = X0[1:8]
        rho = X0[9]        
        mu_Q = X0[10] * hc
        

        P0 = -Omega(phi, Phi1, Phi2, mu_B/hc, mu_Q/hc, T/hc, θ, rho, Gv, ints)
        X_get[i, 2] = P0
        X_get[i, 3:end] = X0
    end

    # 在 θ = π 处：eta 需要变号
    X0[4:6] .*= -1
    println("Calculating segment 2...")
    # 段 2: θ ∈ [π, 3π)
    offset2 = length(θ1)
    for (j, θ) in enumerate(θ2)
        X0 = Tmu(X0, T/hc, θ, mu_B/hc, Gv, ints)
        phi = X0[1:6]
        Phi1 = X0[7]
        Phi2 = X0[8]
        mu_Q = X0[9] * hc
        rho = X0[10]
        P0 = -Omega(phi, Phi1, Phi2, mu_B/hc, mu_Q/hc, T/hc, θ, rho, Gv, ints)
        X_get[offset2 + j, 2] = P0
        X_get[offset2 + j, 3:end] = X0
    end


    # 在 θ = 3π 处：eta 再次变号
    X0[4:6] .*= -1
    println("Calculating segment 3...")
    # 段 3: θ ∈ [3π, 4π]    
    offset3 = length(θ1) + length(θ2)
    for (k, θ) in enumerate(θ3)
        X0 = Tmu(X0, T/hc, θ, mu_B/hc, Gv, ints)
        phi = X0[1:6]
        Phi1 = X0[7]
        Phi2 = X0[8]
        mu_Q = X0[9] * hc
        rho = X0[10]
        P0 = -Omega(phi, Phi1, Phi2, mu_B/hc, mu_Q/hc, T/hc, θ, rho, Gv, ints)
        X_get[offset3 + k, 2] = P0
        X_get[offset3 + k, 3:end] = X0
    end
        # 写出
    data = DataFrame(
        X_get,
        [:theta, :P0, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2, :mu_B, :mu_Q]
    )
    out_file = "../../data/axion_red_vec/vac_sol_Gv=$Gv.csv"
    CSV.write(out_file, data)
    println("Data written to ", out_file)

end

