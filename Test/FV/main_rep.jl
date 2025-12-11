using Revise
includet("../../src/Finite_Vol/Rep.jl")
using CSV 
using DataFrames
using Dierckx


function parametrize_deformation(R, δ;para=2.0,scale=1.0)
    """
    δ: 变形参数 (0 ≤ δ < ∞)
    para: 调节变形幅度的参数
    scale: +1 压扁 -1 拉长
    - δ = 0: 球形 (a=b=c=R)
    - δ > 0 (且 scale=1.0): 扁平椭球(a=b > c)
    - δ > 0 (且 scale=-1.0): 拉长椭球(a=b < c)
    - 表面积单调递增
    """
    V = (4/3)*π*R^3
    
    # 基于 β₂ 的简化
    β₂ = tanh(δ)  # 保证 β₂ < 1
  #@  para = 1.8
    a = R * (1 + para * β₂)^(scale * 2/3)
    b = a
    c = V / ((4/3)*π*a*b)
    
    return a, b, c
end



function main_rep(; mu_B::Float64 = 0.0)
    # mu_BC = 391*3 MEV
    firstline_path = "../../data/FV/1st/first_R30_e00.csv"
    Ts = 10.0:2.0:300.0   # 单位：MeV
    R = 30.0  # fm
    e = 0.0  # 变形参数
    a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)
    ints = get_nodes_el(200, a, b, c, modes="D")

    # ===== 两套典型初值 =====
    X_CONF   = [-1.9, -1.9, -2.2, 0.0038, 0.0038]
    X_DECONF = [-0.5, -0.5, -1.8, 0.80,   0.80  ]
    
    # 读取一阶线表
    df_in = CSV.read(firstline_path, DataFrame)
    df_unique = unique(df_in, :mu_star)
    sort!(df_unique, :mu_star)

    mu_unique = df_unique.mu_star
    T_unique = df_unique.T
    
    # 创建 μ → T 的三次样条插值
    T_of_mu = Spline1D(mu_unique, T_unique, k=3, bc="error")
    
    # 获取插值有效范围
    mu_min, mu_max = extrema(mu_unique)
    T_CEP = maximum(T_unique)  # 临界点温度
    
    # 判断当前 mu_B 是否在一阶线范围内
    in_firstline_range = (mu_min <= mu_B <= mu_max)
    
    # 精确计算对应的临界温度
    T_critical = if in_firstline_range
        T_of_mu(mu_B)
    else
        -Inf  # 不在范围内，设为无穷小
    end
    
    println("mu_B = $mu_B MeV")
    if in_firstline_range
        println("插值得到的临界温度: T_crit = $(round(T_critical, digits=2)) MeV")
    else
        println("mu_B 超出一阶线范围 [$mu_min, $mu_max] MeV，无一阶相变")
    end

    # ===== 判断第一个点 (mu_B, Ts[1]) 的初始相态 =====
    T_start = Ts[1]
    
    # 相变线左下是禁闭相，右上是去禁闭相
    # 左下: T < T_critical(μ) 或 μ < μ_critical(T)
    # 右上: T > T_critical(μ) 或 μ > μ_critical(T)
    initial_confined = if in_firstline_range
        # 在一阶线范围内：直接比较温度
        T_start < T_critical
    else
        # 超出一阶线范围
        if mu_B < mu_min
            # mu 太小：cross over
            true
        elseif mu_B > mu_max
            # mu 太大：判断温度
            false  # 低于临界点温度为禁闭相
        else
            # 不应到达这里
            true
        end
    end
    X00 = [-1.8, -1.8, -2.2, 0.01, 0.01]  # phi_u, phi_d, phi_s, Phi1, Phi2
    X00 = Tmu(50.0/hc, 0.0/hc, X00, ints)
    phi = X00[1:3]
    Phi1 = X00[4]
    Phi2 = X00[5]
    P0 = -Omega(phi, Phi1, Phi2, 50.0/hc, 0.0/hc, ints)
    println("P0=$P0")




    # 根据判断选择初始猜测
    X0 = if initial_confined
        println("初始点在相变线左下（禁闭相区域）")
        copy(X_CONF)
    else
        println("初始点在相变线右上（去禁闭相区域）")
        copy(X_DECONF)
    end

    data = zeros(length(Ts), 7)
    NewX = similar(X0)
    k_1st = 0
    sol = zeros(length(Ts), 6)
    sol[:, 1] .= Ts
    
    for (i, T) in enumerate(Ts)
        println("T = $T")

        # ===== 精确判据：穿越相变线时切换初值 =====
        if in_firstline_range
            # 检测是否穿越相变线
            if initial_confined && T >= T_critical
                # 从禁闭相进入去禁闭相
                k_1st += 1
                if k_1st == 1
                    X0 .= X_DECONF
                    println("  [T ≥ T_crit，穿越相变线，切换到去禁闭相]")
                end
            elseif !initial_confined && T < T_critical
                # 从去禁闭相进入禁闭相
                k_1st += 1
                if k_1st == 1
                    X0 .= X_CONF
                    println("  [T < T_crit，穿越相变线，切换到禁闭相]")
                end
            end
        end

        # 牛顿迭代
        NewX = Tmu(T/hc, mu_B/hc, X0, ints)
        sol[i, 2:end] = NewX
        
        # 微扰计算
        data[i, :] = Ther_Rep(NewX*1.0001, T/hc, mu_B/hc, ints, P0)
        # 延拓
        X0 .= NewX
    end

    df = DataFrame(data, [:T, :mu_B, :P, :E, :CV, :CP, :cs2])
    df_sol = DataFrame(sol, [:T, :phiu, :phid, :phis, :Phi, :Phi_bar])
    CSV.write("../../data/FV/rep_mu=R$(R)e$(e)_$(mu_B).csv", df)
    CSV.write("../../data/FV/rep_sol_mu=R$(R)e$(e)_$(mu_B).csv", df_sol)
    println("Done!")
end


function main_rep_mu(; T::Float64 = 100.0)
    """
    固定温度 T，化学势 mu_B 从小到大扫描
    自动识别并跳过一阶相变线
    """
  

    paths = [
        "../../data/FV/1st/first_R30_e00.csv",
        "../../data/FV/1st/first_R30_e03.csv",
        "../../data/FV/1st/first_R30_e07.csv",
        "../../data/FV/1st/first_R30_e10.csv",
    ]
    Res = [
        (30.0, 0.0),
        (30.0, 0.3),
        (30.0, 0.7),
        (30.0, 1.0),
    ]

   """ 
    paths = [
        "../../data/FV/1st/first_R100_e00.csv",
        "../../data/FV/1st/first_R100_e03.csv",
        "../../data/FV/1st/first_R100_e07.csv",
        "../../data/FV/1st/first_R100_e10.csv",
    ]
    Res = [
        (100.0, 0.0),
        (100.0, 0.3),
        (100.0, 0.7),
        (100.0, 1.0),
    ]
    """ 

    for (path, (R, e)) in zip(paths, Res)
        println("Processing firstline file: $path")
        calc_mu(T, path, R, e)
    end


end

function find_sol_0(T::Float64, R::Float64, e::Float64)
    df_in = "../../data/FV/1st/D_V_R=$(R).csv"
    df = CSV.read(df_in, DataFrame)
    row = filter(row -> isapprox(row.T, T; atol=2) && isapprox(row.e, e; atol=1e-3), eachrow(df))
    if nrow(row) == 0
        error("No solution found for T=$T, R=$R, e=$e")
    end
    phi_u = row.phi_u[1]
    phi_d = row.phi_d[1]
    phi_s = row.phi_s[1]
    Phi1  = row.Phi1[1]
    Phi2  = row.Phi2[1]
    return [phi_u, phi_d, phi_s, Phi1, Phi2]
end



function calc_mu(T, firstline_path, R, e)
        
    mu_Bs = 10.0:6.0:1200.0   # 单位：MeV
    a, b, c = parametrize_deformation(R, e; para=3.0, scale=-1.0)
    ints = get_nodes_el(200, a, b, c, modes="D")

    # 两套典型初值
    X_CONF   = [-1.9, -1.9, -2.2, 0.0038, 0.0038]
    X_DECONF = [-0.5, -0.5, -1.8, 0.80,   0.80  ]

    # 读取并处理一阶线数据 (构造 T -> mu 插值)
    df_in = CSV.read(firstline_path, DataFrame)
    df_unique = unique(df_in, :T)
    sort!(df_unique, :T)
    T_unique = df_unique.T
    mu_unique = df_unique.mu_star

    # 如果数据不足或不单调请根据实际情况处理
    mu_of_T = Spline1D(T_unique, mu_unique, k=1, bc="error")
    T_min, T_max = extrema(T_unique)

    in_firstline_range = (T_min <= T <= T_max)
    mu_critical = in_firstline_range ? mu_of_T(T) : Inf

    println("T = $T MeV, in_firstline_range = $in_firstline_range, mu_crit = $mu_critical")

    # 参考压强 P0（同 main_rep 的做法）
    X00 = [-1.8, -1.8, -2.2, 0.01, 0.01]
    X00 = Tmu(50.0/hc, 0.0/hc, X00, ints)
    phi = X00[1:3]; Phi1 = X00[4]; Phi2 = X00[5]
    P0 = -Omega(phi, Phi1, Phi2, 50.0/hc, 0.0/hc, ints)

    X0 = find_sol_0(T, R, e)

    # 准备保存两段数据（相变线之前 / 之后）
    Tb = Float64[]; mub = Float64[]; Pb = Float64[]; Eb = Float64[]; cs2b = Float64[]
    Ta = Float64[]; mua = Float64[]; Pa = Float64[]; Ea = Float64[]; cs2a = Float64[]
    CVb = Float64[]; CPb = Float64[]; CVa = Float64[]; CPa = Float64[]


    crossed = false

    for mu_B in mu_Bs
        println("mu_B = $mu_B")

        # 如果在一阶线范围并且遇到交界，跳过相变（切换初值）
        if in_firstline_range && (!crossed) && (mu_B >= mu_critical)
            crossed = true
            println("  crossing first-order line at mu ≈ $mu_critical, switch initial guess to deconfined")
            X0 .= X_DECONF
        end

        # 求解并微扰
        NewX = Tmu(T/hc, mu_B/hc, X0, ints)
        vals = Ther_Rep(NewX * 1.0001, T/hc, mu_B/hc, ints, P0)
        # 约定 Ther_Rep 返回向量: [T*197.33, mu*197.33, P, E,  CV, CP， v2,]
        P, E, CV, CP, cs2 = vals[3], vals[4], vals[5], vals[6], vals[7]

        # 根据是否已越过相变线存入不同数组
        if !crossed
            push!(Tb, T); push!(mub, mu_B); push!(Pb, P); push!(Eb, E); push!(cs2b, cs2); push!(CVb, CV); push!(CPb, CP)
        else
            push!(Ta, T); push!(mua, mu_B); push!(Pa, P); push!(Ea, E); push!(cs2a, cs2); push!(CVa, CV); push!(CPa, CP)
        end

        # 延拓用于下一步
        X0 .= NewX
    end

    # 写文件：若未穿越则只有 before 文件（或写为 single）
    if crossed
        df_before = DataFrame(T = Tb, mu_B = mub, P = Pb, E = Eb, cs2 = cs2b, CV = CVb, CP = CPb)
        df_after  = DataFrame(T = Ta, mu_B = mua, P = Pa, E = Ea, cs2 = cs2a, CV = CVa, CP = CPa)
        CSV.write("../../data/FV/cs/rep_T=$(round(T,digits=3))_R=$(R)_e=$(e)_before.csv", df_before)
        CSV.write("../../data/FV/cs/rep_T=$(round(T,digits=3))_R=$(R)_e=$(e)_after.csv",  df_after)
        println("Saved two files: before & after crossing.")
    else
        df_all = DataFrame(T = Tb, mu_B = mub, P = Pb, E = Eb, cs2 = cs2b, CV = CVb, CP = CPb)
        CSV.write("../../data/FV/cs/rep_T=$(round(T,digits=3))_R=$(R)_e=$(e)_all.csv", df_all)
        println("No crossing found, saved single file.")
    end
end