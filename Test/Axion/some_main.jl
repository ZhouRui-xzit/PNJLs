include("../../src/axion/pnjl_axion.jl")  # ← 改为 axion 版本的实现文件
using DataFrames
using CSV
using BenchmarkTools



function main_Tmu()
    Ts = 300:-1:100.0
    mu_B = 0.0:2.0*3:300.0*3
    # X0: [sigma_u, sigma_d, sigma_s, eta_u, eta_d, eta_s, Phi1, Phi2]
    X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1]
    lens = length(Ts) * length(mu_B)
    data = zeros(lens, 7)  # T, mu_B, (phi_u→sigma_u), (phi_d→sigma_d), (phi_s→sigma_s), Phi1, Phi2
    for (i, MU) in enumerate(mu_B)
        println("mu_B = $MU MEV")
        X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1]  # 重置初始猜测
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            # axion 接口：Tmu(T_fm^-1, muB_fm^-1, X0, theta)
            X0 = Tmu(T/hc, MU/hc, X0, THETA)
            # 输出列名保持不变：phi_* 列填入 sigma_*
            data[idx, :] = [T, MU, X0[1], X0[2], X0[3], X0[7], X0[8]]
        end
    end
    df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/axion/T_mu_B_scan.dat", df)
end

function main_Trho()
    T1s = 131.03:-0.1:130.03
    T2s = 130:-2:10.0
    Ts = vcat(T1s, T2s)
    rho1s = 3.00:-0.01:0.01
    rho2s = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    rho2s = reverse(rho2s)
    rhos = vcat(rho1s, rho2s)
    # X0: [sigma(3), eta(3), Phi1, Phi2, muB_div3, aux1, aux2]
    X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1, 320/hc, 320/hc, 320/hc]
    lens = length(Ts) * length(rhos)
    data = zeros(lens, 4)  # T, rho_B, mu, P
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1, 320/hc, 320/hc, 320/hc]  # 重置初始猜测
        for (j, rho_B) in enumerate(rhos)
            idx = (i-1)*length(rhos) + j
            # axion 接口：Trho(T_fm^-1, rho_B, X0, theta)
            X0 = Trho(T/hc, rho_B, X0, THETA)
            # 仍沿用原脚本习惯：mu = (mu_B/3) 的 MeV 值
            mu = X0[9] * hc    # X0[9] 是 mu_B/3 (fm^-1)
            # 取出序参量
            phi6 = X0[1:6]     # [sigma(3), eta(3)]
            Phi1 = X0[7]
            Phi2 = X0[8]
            # 压强：保持你原来的“mu/hc”传参方式（最小改动）
            P = -Omega(phi6, Phi1, Phi2, T/hc, mu/hc, THETA)
            data[idx, :] = [T, rho_B, mu, P]
        end
    end
    df = DataFrame(data, [:T, :rho_B, :mu, :P])
    CSV.write("../../data/axion/T_rho_B_scan.dat", df)
end

function main_mu(mu_B)
    Ts = 300:-2:1.0

    thetas = [0.0]
    # X0: 8 维
    X0 = [0.01, 0.01, 0.40, 0.01, 0.01, 0.01, 0.8, 0.8]
    data = zeros(length(Ts), 7)
    for theta in thetas
        for (i, T) in enumerate(Ts)
            println("T = $T, mu_B = $mu_B")
            X0 = Tmu(T/hc, mu_B/hc, X0, theta)
            data[i, :] = [T, mu_B, X0[1], X0[2], X0[3], X0[7], X0[8]]
        end
        df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
        CSV.write("../../data/Axion/mu_B=$mu_B theta=$(round(theta, digits=2)).dat", df)
    end

end

function main_rho(T)
    rhos = 0.1:0.01:3.0
    #rhos2 = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    #rhos2 = reverse(rhos2)
    #rhos = vcat(rhos1, rhos2)

    theta = 0.0

    # X0: 11 维
    X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1, 320/hc, 320/hc, 320/hc]
    #X0 = [0.09441, 0.09437, 2.21244, 1.74724, 1.74724, 0.0683, 0.00154, 0.00183, 1.84653, 1.84653, 1.84653]
    data = zeros(length(rhos), length(X0) + 2)
    for (i, rho_B) in enumerate(rhos)
        println("T = $T, rho_B = $rho_B")
        X0 = Trho(T/hc, rho_B, X0, theta)

        data[i, 1] = T
        data[i, 2] = rho_B
        data[i, 3:end] = X0
    end
    df = DataFrame(data, [:T, :rho_B, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2, :mu_u, :mu_d, :mu_s])
    CSV.write("../../data/axion/T=$T.dat", df)
end

function single_mu()
    T = 300.0
    mu_B = 0.0
    theta = 0.0
    X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1]
    _ = Tmu(T/hc, mu_B/hc, X0, theta)
end

function single_rho()
    T = 10.0
    rho_B = 3.0
    theta = 0.0
    X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1, 320/hc, 320/hc, 320/hc]
    _ = Trho(T/hc, rho_B, X0, theta)
end


function main_theta_mu()
    # 扫描参数
    dθ   = 0.01
    ϵ    = 1e-6  # 避免端点重合
    θ1   = range(0.0,     step=dθ, stop=pi - ϵ)
    θ2   = range(pi,      step=dθ, stop=3pi - ϵ)
    θ3   = range(3pi,     step=dθ, stop=4pi)
    lens = length(θ1) + length(θ2) + length(θ3)

    # 物理参数（外部用 MeV；内部调用时会除以 hc）
    mu_B = 0.0
    T    = 5.0

    # 初值（axion 接口：X0 = [sigma(3), eta(3), Phi1, Phi2]）
    X0   = [1.8, 1.8, 2.0,   0.1, 0.1, 0.1,   0.1, 0.1]

    # 结果容器：列 = [:T, :theta, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2]
    X_get = zeros(lens, length(X0) + 3)
    X_get[:, 1] .= T
    X_get[:, 2] .= vcat(collect(θ1), collect(θ2), collect(θ3))

    # 输出文件
    output_file = "../../data/axion/axion_thetas.dat"

    # 段 1: θ ∈ [0, π)
    println("计算 θ ∈ [0, π) ...")
    for (i, θ) in enumerate(θ1)
        X0 = Tmu(T/hc, mu_B/hc, X0, θ)
        P0 = - Omega(X0[1:6], X0[7], X0[8], T/hc, mu_B/hc, θ)
        X_get[i, 3] = P0 
        X_get[i, 4:end] = X0
    end

    # 在 θ = π 处：eta 需要变号
    X0[4:6] .*= -1

    # 段 2: θ ∈ [π, 3π)
    offset2 = length(θ1)
    println("计算 θ ∈ [π, 3π) ...")
    for (j, θ) in enumerate(θ2)
        X0 = Tmu(T/hc, mu_B/hc, X0, θ)
        P0 = - Omega(X0[1:6], X0[7], X0[8], T/hc, mu_B/hc, θ)
        X_get[offset2 + j, 3] = P0
        X_get[offset2 + j, 4:end] = X0
    end

    # 在 θ = 3π 处：eta 再次变号
    X0[4:6] .*= -1

    # 段 3: θ ∈ [3π, 4π]
    println("计算 θ ∈ [3π, 4π] ...")
    offset3 = length(θ1) + length(θ2)
    for (k, θ) in enumerate(θ3)
        X0 = Tmu(T/hc, mu_B/hc, X0, θ)
        P0 = - Omega(X0[1:6], X0[7], X0[8], T/hc, mu_B/hc, θ)
        X_get[offset3 + k, 3] = P0
        X_get[offset3 + k, 4:end] = X0
    end

    # 写出
    data = DataFrame(
        X_get,
        [:T, :theta, :P0, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2]
    )
    CSV.write(output_file, data)
end

function main_theta_rho()
    # 扫描参数
    dθ   = 0.01
    ϵ    = 1e-6  # 避免端点重合
    θ1   = range(0.0,     step=dθ, stop=pi - ϵ)
    θ2   = range(pi,      step=dθ, stop=3pi - ϵ)
    θ3   = range(3pi,     step=dθ, stop=4pi)
    lens = length(θ1) + length(θ2) + length(θ3)

    # 物理参数（外部用 MeV；内部调用时会除以 hc）
    mu_B = 0.0
    T    = 5.0
    rho  = 0.1
    # 初值（axion 接口：X0 = [sigma(3), eta(3), Phi1, Phi2]）
    X0   = [1.8, 1.8, 2.0,   0.1, 0.1, 0.1,   0.1, 0.1, 300/hc, 300/hc, 300/hc]

    # 结果容器：列 = [:rho, :theta, :P, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2]
    X_get = zeros(lens, length(X0) + 3)
    X_get[:, 1] .= rho
    X_get[:, 2] .= vcat(collect(θ1), collect(θ2), collect(θ3))

    # 输出文件
    output_file = "../../data/axion/axion_thetas_rho.dat"

    # 段 1: θ ∈ [0, π)
    println("计算 θ ∈ [0, π) ...")
    for (i, θ) in enumerate(θ1)
        println("theta = $θ")
        X0 = Trho(T/hc, rho, X0, θ)
        muB = X0[9] * 3
        P0 = - Omega(X0[1:6], X0[7], X0[8], T/hc, muB/hc, θ)
        X_get[i, 3] = P0 
        X_get[i, 4:end] = X0
    end

    # 在 θ = π 处：eta 需要变号
    X0[4:6] .*= -1

    # 段 2: θ ∈ [π, 3π)
    offset2 = length(θ1)
    println("计算 θ ∈ [π, 3π) ...")
    for (j, θ) in enumerate(θ2)
        X0 = Trho(T/hc, rho, X0, θ)
        muB = X0[9] * 3
        P0 = - Omega(X0[1:6], X0[7], X0[8], T/hc, muB/hc, θ)
        X_get[offset2 + j, 3] = P0
        X_get[offset2 + j, 4:end] = X0
    end

    # 在 θ = 3π 处：eta 再次变号
    X0[4:6] .*= -1

    # 段 3: θ ∈ [3π, 4π]
    println("计算 θ ∈ [3π, 4π] ...")
    offset3 = length(θ1) + length(θ2)
    for (k, θ) in enumerate(θ3)
        X0 = Trho(T/hc, rho, X0, θ)
        muB = X0[9] * 3
        P0 = - Omega(X0[1:6], X0[7], X0[8], T/hc, muB/hc, θ)
        X_get[offset3 + k, 3] = P0
        X_get[offset3 + k, 4:end] = X0
    end

    # 写出
    data = DataFrame(
        X_get,
        [:rho, :theta, :P0, :sigma_u, :sigma_d, :sigma_s, :eta_u, :eta_d, :eta_s, :Phi1, :Phi2,:mu_u, :mu_d, :mu_s]
    )
    CSV.write(output_file, data)
end






if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
