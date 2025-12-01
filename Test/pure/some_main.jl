include("../../src/pure/Pnjl_pure.jl")
using DataFrames
using CSV
using BenchmarkTools



function main_Tmu()
    Ts = 300:-1:100.0
    mu_B = 0.0:2.0*3:300.0*3
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    lens = length(Ts) * length(mu_B)
    data = zeros(lens, 7)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, MU) in enumerate(mu_B)
        println("mu_B = $MU MEV")
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # 重置初始猜测
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            
            X0 = Tmu(T/hc, MU/hc, X0)
            data[idx, :] = [T, MU, X0...]
        end
    end
    df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/pure/T_mu_B_scan.dat", df)
end

function main_Trho()

    ints1 = get_nodes(200, nodes2=200)
    ints2 = get_nodes(256, nodes2=500)

    T1s = range(131.015625, 131.01, length=10)
    T2s = range(131.00, 130.00, length=10)
    T3s = range(130.0, 10.0, length=50)
    Ts = unique!(vcat(collect(T1s), collect(T2s), collect(T3s)))

    rho1s = 3.00:-0.01:0.01
    rho2s = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    rho2s = reverse(rho2s)
    rhos = vcat(rho1s, rho2s)
    X00 = [-1.8,-1.8, -2.2, 0.01,0.01, 1000/hc]  # 重置初始猜测
    X0 = similar(X00)
    lens = length(Ts) * length(rhos)
    data = zeros(lens, 8)  # T, rho_B, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        if T >= 30.0
            ints = ints1
        else
            ints = ints2
        end
        X0 = X00
        for (j, rho_B) in enumerate(rhos)
            idx = (i-1)*length(rhos) + j
            
            X0 = Trho(T/hc, rho_B, X0, ints)
            mu = X0[6] * hc 
            
            data[idx, :] = [T, rho_B, mu, X0[1:5]...]
            if j==1
                X00 = X0
            end
        end
    end
    df = DataFrame(data, [:T, :rho_B, :mu, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/pure/T_rho_B_scan1123.dat", df)
end





function main_mu(mu_B)
    Ts = 300:-1:0.1
    
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    data = zeros(length(Ts), 7)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, T) in enumerate(Ts)
        println("T = $T, mu_B = $mu_B")
        X0 = Tmu(T/hc, mu_B/hc, X0)
        data[i, :] = [T, mu_B, X0...]
    end
    df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/pure/mu_B=$mu_B.dat", df)
end

function main_rho(T)
    rhos1 = 3.00:-0.01:0.01
    rhos2 = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    rhos2 = reverse(rhos2)
    rhos = vcat(rhos1, rhos2)
    
    X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 320/hc, 320/hc, 320/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    data = zeros(length(rhos), 3) 
    for (i, rho_B) in enumerate(rhos)
        println("T = $T, rho_B = $rho_B")
        X0 = Trho(T/hc, rho_B, X0)
        mu = X0[6] * hc 
        data[i, :] = [T, rho_B, mu]
    end
    df = DataFrame(data, [:T, :rho_B, :mu_B])
    CSV.write("../../data/pure/T=$T.dat", df)
end

function single_mu()
    T = 300.0
    mu_B = 0.0
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    X_sol = Tmu(T/hc, mu_B/hc, X0)
    #println("T = $T, mu_B = $mu_B")
    #println("Solution: phi_u=$(X_sol[1]), phi_d=$(X_sol[2]), phi_s=$(X_sol[3]), Phi1=$(X_sol[4]), Phi2=$(X_sol[5])")
end

function single_rho()
    T = 10.0
    rho_B = 3.0
    X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 320/hc, 320/hc, 320/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    X_sol = Trho(T/hc, rho_B, X0)
    mu = X_sol[6] * hc 
    #println("T = $T, rho_B = $rho_B")
    #println("Solution: phi_u=$(X_sol[1]), phi_d=$(X_sol[2]), phi_s=$(X_sol[3]), Phi1=$(X_sol[4]), Phi2=$(X_sol[5]), mu_B=$mu")
end

function main()
    while true
        println("\n===== PNJL 模型计算程序 =====")
        println("请选择要运行的程序：")
        println("  1) 定 μ_B 扫描 T (调用 main_mu)")
        println("  2) 定 T 扫描 ρ_B (调用 main_rho)")
        println("  3) 扫描 T-μ_B 二维参数空间 (调用 main_Tmu)")
        println("  4) 扫描 T-ρ_B 二维参数空间 (调用 main_Trho)")
        println("  5) 计算单个 T-μ_B 点 (调用 single_mu)")
        println("  6) 计算单个 T-ρ_B 点 (调用 single_rho)")
        println("  7) 性能测试 (使用 BenchmarkTools)")
        println("  q) 退出")
        print("输入选项并回车: "); flush(stdout)
        choice = try
            chomp(readline(stdin))
        catch
            return
        end
        c = lowercase(choice)
        if c in ("q", "quit", "exit")
            println("已退出。")
            return
        elseif c == "1"
            print("请输入 μ_B (单位：MeV): "); flush(stdout)
            s = try
                chomp(readline(stdin))
            catch e
                ""
            end
            try
                muB = parse(Float64, s)
                main_mu(muB)
                println("\n计算完成! 结果已保存到: ../../data/pure/mu_B=$muB.dat")
            catch e
                println("计算出错: $e")
            end
        elseif c == "2"
            print("请输入 T (单位：MeV): "); flush(stdout)
            s = try chomp(readline(stdin)) catch e ""; end
            try
                T = parse(Float64, s)
                main_rho(T)
                println("\n计算完成! 结果已保存到: ../../data/pure/T=$T.dat")
            catch e
                println("计算出错: $e")
            end
        elseif c == "3"
            println("执行 T-μ_B 二维参数扫描...")
            try
                main_Tmu()
                println("\n计算完成! 结果已保存到: ../../data/pure/T_mu_B_scan.dat")
            catch e
                println("计算出错: $e")
            end
        elseif c == "4"
            println("执行 T-ρ_B 二维参数扫描...")
            try
                main_Trho()
                println("\n计算完成! 结果已保存到: ../../data/pure/T_rho_B_scan.dat")
            catch e
                println("计算出错: $e")
            end
        elseif c == "5"
            println("计算单个 T-μ_B 点...")
            print("请输入 T (单位：MeV) [默认300.0]: "); flush(stdout)
            s_T = try chomp(readline(stdin)) catch e ""; end
            print("请输入 μ_B (单位：MeV) [默认0.0]: "); flush(stdout)
            s_mu = try chomp(readline(stdin)) catch e ""; end
            
            try
                T = isempty(s_T) ? 300.0 : parse(Float64, s_T)
                mu_B = isempty(s_mu) ? 0.0 : parse(Float64, s_mu)
                
                X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]
                X_sol = Tmu(T/hc, mu_B/hc, X0)
                
                println("\n计算结果:")
                println("  T = $T MeV, μ_B = $mu_B MeV")
                println("  φ_u = $(X_sol[1]), φ_d = $(X_sol[2]), φ_s = $(X_sol[3])")
                println("  Φ_1 = $(X_sol[4]), Φ_2 = $(X_sol[5])")
            catch e
                println("计算出错: $e")
            end
        elseif c == "6"
            println("计算单个 T-ρ_B 点...")
            print("请输入 T (单位：MeV) [默认10.0]: "); flush(stdout)
            s_T = try chomp(readline(stdin)) catch e ""; end
            print("请输入 ρ_B [默认3.0]: "); flush(stdout)
            s_rho = try chomp(readline(stdin)) catch e ""; end
            
            try
                T = isempty(s_T) ? 10.0 : parse(Float64, s_T)
                rho_B = isempty(s_rho) ? 3.0 : parse(Float64, s_rho)
                
                X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 320/hc, 320/hc, 320/hc]
                X_sol = Trho(T/hc, rho_B, X0)
                mu = X_sol[6] * hc
                
                println("\n计算结果:")
                println("  T = $T MeV, ρ_B = $rho_B")
                println("  φ_u = $(X_sol[1]), φ_d = $(X_sol[2]), φ_s = $(X_sol[3])")
                println("  Φ_1 = $(X_sol[4]), Φ_2 = $(X_sol[5])")
                println("  μ_B = $mu MeV")
            catch e
                println("计算出错: $e")
            end
        elseif c == "7"
            println("执行性能测试...")
            println("\n测试 single_mu 函数:")
            @btime single_mu()
            println("\n测试 single_rho 函数:")
            @btime single_rho()
        else
            println("无效选项: $choice")
        end
        println()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_mu(0.0)
end