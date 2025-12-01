using Revise
using CSV, DataFrames
using Plots

includet("../../src/axion_red_vec/pnjl_vec.jl")

function generate_rhos2(n::Int)
    rhos2 = Float64[]
    for e in -n:-2
        push!(rhos2, 1.0 * 10.0^e)
        push!(rhos2, 5.0 * 10.0^e)
    end
    #push!(rhos2, 0.05)  # 确保包含终点
    return sort(rhos2)
end


function data_Gv(;theta=0.0)
    # theta = 0.0 P0== 25.38212212875826
    sol_path = "../../data/axion_red_vec/vac_sol_Gv=0.5.csv"
    sol = CSV.read(sol_path, DataFrame)
    # 单行找到最近索引
    idx = argmin(abs.(sol[:, 1] .- theta))
    println("Using theta = ", sol[idx, 1])
    theta = sol[idx, 1]
    X0 = Vector(sol[idx, 3:end])
    X0[end-1:end] .*= 1/hc  # to fm^-4
    P0 = sol[idx, 2]

    T = 5.0
    rhos1 = 0.1:0.1:25.0
    rhos2 = generate_rhos2(11)
    rhos = sort(vcat(collect(rhos1), rhos2))
    


    Rv = 0.5
    Gv = Rv * G_f

   
    ints = get_nodes(800)

    data = zeros(length(rhos), 14)
    EOS = zeros(length(rhos), 3)
    data[:, 1] = rhos
    EOS[:, 1] = rhos
    for (i, rho) in enumerate(rhos)
        #println("Calculating rho = ", rho)
        X0 = Trho(X0, T/hc ,rho, theta, Gv, ints)
        data[i, 2:end] = X0
      
        
        EOS[i, 2:3] = get_EOS(X0, T/hc, theta, Gv, P0, ints)
        if EOS[i, 2] < 0 
            println("Negative pressure at rho = ", rho)
        end
    end
    EOS[:, end-1:end] .*= hc # to Mev/fm3
    data[:, end-1:end] .*= hc # to Mev/fm3
    df = DataFrame(data, [:rho, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2, :rho_u, :rho_d, :rho_s, :mu_B, :mu_Q])
    df2 = DataFrame(EOS, [:rho, :P, :E])
    
    outpath1 = "../../data/axion_red_vec/Rv=$(Rv)_theta=$theta.csv"
    outpath2 = "../../data/axion_red_vec/EOS_Rv=$(Rv)_theta=$theta.csv"
    CSV.write(outpath1, df)
    CSV.write(outpath2, df2)

end



function find_P0_mu()
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
        mu_u = 1/3*mu_B + 2/3 * mu_Q
        mu_d = 1/3*mu_B - 1/3 * mu_Q
        mu_s = 1/3*mu_B - 1/3 * mu_Q
        mu_e = - mu_Q
        mu_mu = - mu_Q
        mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
        P0 = -Omega(orders, mus/hc, T/hc, θ, rho, Gv, ints)
        X_get[i, 2] = P0
        X_get[i, 3:end] = X0
    end

    # 在 θ = π 处：eta 需要变号
    X0[4:6] .*= -1
    println("Calculating segment 2...")
    # 段 2: θ ∈ [π, 3π)
    offset2 = length(θ1)
    for (j, θ) in enumerate(θ2)
        X0 = Tmu(X0, mu_B/hc, T/hc, θ, Gv, ints)
        orders = X0[1:8]
        rho = X0[9]        
        mu_Q = X0[10] * hc
        mu_u = 1/3*mu_B + 2/3 * mu_Q
        mu_d = 1/3*mu_B - 1/3 * mu_Q
        mu_s = 1/3*mu_B - 1/3 * mu_Q
        mu_e = - mu_Q
        mu_mu = - mu_Q
        mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
        P0 = -Omega(orders, mus/hc, T/hc, θ, rho, Gv, ints)
        X_get[offset2 + j, 2] = P0
        X_get[offset2 + j, 3:end] = X0
    end


    # 在 θ = 3π 处：eta 再次变号
    X0[4:6] .*= -1
    println("Calculating segment 3...")
    # 段 3: θ ∈ [3π, 4π]    
    offset3 = length(θ1) + length(θ2)
    for (k, θ) in enumerate(θ3)
        X0 = Tmu(X0, mu_B/hc, T/hc, θ, Gv, ints)
        orders = X0[1:8]
        rho = X0[9]        
        mu_Q = X0[10] * hc
        mu_u = 1/3*mu_B + 2/3 * mu_Q
        mu_d = 1/3*mu_B - 1/3 * mu_Q
        mu_s = 1/3*mu_B - 1/3 * mu_Q
        mu_e = - mu_Q
        mu_mu = - mu_Q
        mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
        P0 = -Omega(orders, mus/hc, T/hc, θ, rho, Gv, ints)
        X_get[offset3 + k, 2] = P0
        X_get[offset3 + k, 3:end] = X0
    end
    X_get[:, end] .*= hc

    data = DataFrame(
        X_get,
        [:theta, :P0, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2, :rho, :mu_Q]
    )
    out_file = "../../data/axion_red_vec/vac_sol_Gv=$Rv.csv"
    CSV.write(out_file, data)
    println("Data written to ", out_file)

end



function find_P0_rho()
    T = 5.0
    rho = 1e-12
    dθ   = 0.01
    ϵ    = 1e-6  # 避免端点重合
    Rv = 0.500
    Gv = Rv * G_f
    θ1   = range(0.0,     step=dθ, stop=pi - ϵ)
    θ2   = range(pi,      step=dθ, stop=3pi - ϵ)
    θ3   = range(3pi,     step=dθ, stop=4pi)
    lens = length(θ1) + length(θ2) + length(θ3)


   # mu_B = 1122.78/hc
   # mu_Q = -7.07/hc


    X0 = [1.8, 1.8, 2.2, -0.001, -0.001, -0.001, 0.00383, 0.00383, 0.01, 0.01, 0.01, 1000.0/hc, -10.0/hc]
 
    X_get = zeros(lens, length(X0) + 2)
  
    X_get[:, 1] .= vcat(collect(θ1), collect(θ2), collect(θ3))

    ints = get_nodes(256)

     # 段 1: θ ∈ [0, π)
    println("Calculating segment 1...")
    for (i, θ) in enumerate(θ1)
        X0 = Trho(X0, T/hc, rho, θ, Gv, ints)
        orders = X0[1:8]
        rhos = X0[9:11]
        mu_B = X0[12] * hc
        mu_Q = X0[13] * hc
        mu_u = 1/3*mu_B + 2/3 * mu_Q
        mu_d = 1/3*mu_B - 1/3 * mu_Q
        mu_s = 1/3*mu_B - 1/3 * mu_Q
        mu_e = - mu_Q
        mu_mu = - mu_Q
        mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
        P0 = -Omega(orders, mus/hc, T/hc, θ, rhos, Gv, ints)
        X_get[i, 2] = P0
        X_get[i, 3:end] = X0
    end

    # 在 θ = π 处：eta 需要变号
    X0[4:6] .*= -1
    println("Calculating segment 2...")
    # 段 2: θ ∈ [π, 3π)
    offset2 = length(θ1)
    for (j, θ) in enumerate(θ2)
        X0 = Trho(X0, T/hc, rho, θ, Gv, ints)
        orders = X0[1:8]
        rhos = X0[9:11]
        mu_B = X0[12] * hc
        mu_Q = X0[13] * hc
        mu_u = 1/3*mu_B + 2/3 * mu_Q
        mu_d = 1/3*mu_B - 1/3 * mu_Q
        mu_s = 1/3*mu_B - 1/3 * mu_Q
        mu_e = - mu_Q
        mu_mu = - mu_Q
        mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
        P0 = -Omega(orders, mus/hc, T/hc, θ, rhos, Gv, ints)
        X_get[offset2 + j, 2] = P0
        X_get[offset2 + j, 3:end] = X0
    end


    # 在 θ = 3π 处：eta 再次变号
    X0[4:6] .*= -1
    println("Calculating segment 3...")
    # 段 3: θ ∈ [3π, 4π]    
    offset3 = length(θ1) + length(θ2)
    for (k, θ) in enumerate(θ3)
        X0 = Trho(X0, T/hc, rho, θ, Gv, ints)
        orders = X0[1:8]
        rhos = X0[9:11]
        mu_B = X0[12] * hc
        mu_Q = X0[13] * hc
        mu_u = 1/3*mu_B + 2/3 * mu_Q
        mu_d = 1/3*mu_B - 1/3 * mu_Q
        mu_s = 1/3*mu_B - 1/3 * mu_Q
        mu_e = - mu_Q
        mu_mu = - mu_Q
        mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
        P0 = -Omega(orders, mus/hc, T/hc, θ, rhos, Gv, ints)
        X_get[offset3 + k, 2] = P0
        X_get[offset3 + k, 3:end] = X0
    end
    X_get[:, end-1:end] .*= hc

    data = DataFrame(
        X_get,
        [:theta, :P0, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2, :rho_u, :rho_d, :rho_s, :mu_B, :mu_Q]
    )
    out_file = "../../data/axion_red_vec/vac_sol_Gv=$Rv.csv"
    CSV.write(out_file, data)
    println("Data written to ", out_file)

end

function find_P0_Gv()
    T = 5.0
    mu_B = 1016.58349/hc
    mu_Q = -7.43395E-7/hc
    rho = 1e-11
  
    #1.84321	1.84321	2.22691	4.68e-321	2.925e-321	3.7e-322	1.00018E-7	1.36046E-5	1045.10088	-0.00275
    X0 = [1.8, 1.8, 2.22691, -0.00001, -0.00001, -0.00001, 0.0000383, 0.0000383, mu_B, mu_Q]
    Rvs = 0:0.025:1.0

    theta = 0.0
    lens = length(Rvs)
    X_get = zeros(lens, length(X0) + 2)
    X_get[:, 1] .= Rvs
    ints = get_nodes(800)
    for Rv in Rvs
        println("Calculating Rv = ", Rv)
        Gv = Rv * G_f
        X0 = Trho(X0, T/hc, rho, theta, Gv, ints)
        orders = X0[1:8]
        mu_B = X0[9] * hc
        mu_Q = X0[10] * hc
        mu_u = 1/3*mu_B + 2/3 * mu_Q
        mu_d = 1/3*mu_B - 1/3 * mu_Q
        mu_s = 1/3*mu_B - 1/3 * mu_Q
        mu_e = - mu_Q
        mu_mu = - mu_Q
        mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
        P0 = -Omega(orders, mus/hc, T/hc, theta, rho, Gv, ints)
        X_get[Rv .== Rvs, 2] .= P0
        X_get[Rv .== Rvs, 3:end] = X0
    end
    X_get[:, end-1:end] .*= hc
    data = DataFrame(
        X_get,
        [:Rv, :P0, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2, :mu_B, :mu_Q]
    )
    out_file = "../../data/axion_red_vec/vac_sol_Trho.csv"
    CSV.write(out_file, data)
    println("Data written to ", out_file)
end