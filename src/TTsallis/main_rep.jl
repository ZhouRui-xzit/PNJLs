using Revise
includet("Rep.jl")
using CSV 
using DataFrames
using Dierckx




function main_rep_T(; mu_B::Float64 = 0.0)
    """
    固定温度 T，化学势 mu_B 从小到大扫描
    自动识别并跳过一阶相变线
    """
  

    paths = [
        "1st/first_q_1.001.csv",
        "1st/first_q_1.001.csv",
        "1st/first_q_1.001.csv",

    ]
    Res = [
        1.0001, 
        1.05, 
        1.10
    ]





    for (path, q) in zip(paths, Res)
        println("Processing firstline file: $path")
        calc_T(mu_B, path, q)
        #calc_flu(mu_B, path, q)
    end
end




function calc_T(mu_B::Float64, firstline_path::String, q::Float64)
    """
    固定 mu_B,扫描温度 T
    """
    Ts = 50.0:2.0:500.0   # 单位:MeV
   # Ts2 = 0.1:0.2:1.0
    #Ts = unique(vcat(collect(Ts2), collect(Ts)))
    ints = get_nodes(10; nodes2=64)

    # 两套典型初值
    X_CONF   = [-1.9, -1.9, -2.2, 0.0038, 0.0038]
    X_DECONF = [-0.5, -0.5, -1.8, 0.80,   0.80  ]

    # 读取一阶线数据 (构造 mu -> T 插值)
    df_in = CSV.read(firstline_path, DataFrame)
    df_unique = unique(df_in, :mu_star)
    sort!(df_unique, :mu_star)
    mu_unique = df_unique.mu_star
    T_unique = df_unique.T

    T_of_mu = Spline1D(mu_unique, T_unique, k=1, bc="error")
    mu_min, mu_max = extrema(mu_unique)

    in_firstline_range = (mu_min <= mu_B <= mu_max)
    T_critical = in_firstline_range ? T_of_mu(mu_B) : -Inf

    println("mu_B = $mu_B MeV, in_firstline_range = $in_firstline_range, T_crit = $T_critical")

    # 参考压强
    X00 = [-1.8, -1.8, -2.2, 0.01, 0.01]
    X00 = Tmu(0.1/hc, 0.0/hc, X00, q, ints)
    phi = X00[1:3]; Phi1 = X00[4]; Phi2 = X00[5]
    P0 = -Omega(phi, Phi1, Phi2, 0.1/hc, 0.0/hc, q, ints)
    println("Reference pressure P0 = $P0 MeV/fm^3")
    # 判断起始相态
    T_start = Ts[1]
    initial_confined = in_firstline_range ? (T_start < T_critical) : (mu_B < mu_min)
    X0 = initial_confined ? copy(X_CONF) : copy(X_DECONF)
    println("Initial phase: $(initial_confined ? "confined" : "deconfined")")

    # 准备数据存储
    Tb = Float64[]; mub = Float64[]; Pb = Float64[]; Eb = Float64[]; TAa = Float64[] # TA:= trace anomaly
    Ta = Float64[]; mua = Float64[]; Pa = Float64[]; Ea = Float64[]; TAb = Float64[]
    Sa = Float64[]; rhoa = Float64[]; cs2a = Float64[]; CVa = Float64[];
    Sb = Float64[]; rhob = Float64[]; cs2b = Float64[]; CVb = Float64[];

    crossed = false

    for T in Ts
        println("T = $T")

        # 检测穿越相变线
        if in_firstline_range && (!crossed) && (T >= T_critical)
            crossed = true
            println("  crossing first-order line at T ≈ $T_critical, switching to deconfined phase")
            X0 .= X_DECONF
        end

        # 求解
        NewX = Tmu(T/hc, mu_B/hc, X0, q, ints)
        vals = Ther_Rep(NewX * 1.0001, T/hc, mu_B/hc, P0, q, ints)
        P, E, TA, S, rhoB, Cv, cs2 = vals[3], vals[4], vals[5], vals[6], vals[7], vals[8], vals[9]

        # 存储数据
        if !crossed
            push!(Tb, T); push!(mub, mu_B); push!(Pb, P); push!(Eb, E); push!(TAb, TA);
            push!(Sb, S); push!(rhob, rhoB); push!(cs2b, cs2)
            push!(CVb, Cv)
        else
            push!(Ta, T); push!(mua, mu_B); push!(Pa, P); push!(Ea, E); push!(TAa, TA);
            push!(Sa, S); push!(rhoa, rhoB); push!(cs2a, cs2)
            push!(CVa, Cv)
        end

        X0 .= NewX
    end

    # 写文件
    if crossed
        df_before = DataFrame(T=Tb, mu_B=mub, P=Pb, E=Eb, TA=TAb, S=Sb, rhoB=rhob, cs2=cs2b, CV=CVb)
        df_after  = DataFrame(T=Ta, mu_B=mua, P=Pa, E=Ea, TA=TAa, S=Sa, rhoB=rhoa, cs2=cs2a, CV=CVa)
        CSV.write("cs/rep_muB=$(round(mu_B,digits=3))_q=$(q)_before.csv", df_before)
        CSV.write("cs/rep_muB=$(round(mu_B,digits=3))_q=$(q)_after.csv",  df_after)
        println("Saved two files: before & after crossing.")
    else
        df_all = DataFrame(T=Tb, mu_B=mub, P=Pb, E=Eb, TA=TAb, S=Sb, rhoB=rhob, cs2=cs2b, CV=CVb)
        CSV.write("cs/rep_muB=$(round(mu_B,digits=3))_q=$(q)_all.csv", df_all)
        println("No crossing found, saved single file.")
    end
end


function calc_flu(mu_B::Float64, firstline_path::String, q::Float64)
        """
    固定 mu_B,扫描温度 T
    """
    Ts = 10.0:2.0:300.0   # 单位:MeV
   # Ts2 = 0.1:0.2:1.0
    #Ts = unique(vcat(collect(Ts2), collect(Ts)))
    ints = get_nodes(128, nodes2=128)

    # 两套典型初值
    X_CONF   = [-1.9, -1.9, -2.2, 0.0038, 0.0038]
    X_DECONF = [-0.5, -0.5, -1.8, 0.80,   0.80  ]

    # 读取一阶线数据 (构造 mu -> T 插值)
    df_in = CSV.read(firstline_path, DataFrame)
    df_unique = unique(df_in, :mu_star)
    sort!(df_unique, :mu_star)
    mu_unique = df_unique.mu_star
    T_unique = df_unique.T

    T_of_mu = Spline1D(mu_unique, T_unique, k=1, bc="error")
    mu_min, mu_max = extrema(mu_unique)

    in_firstline_range = (mu_min <= mu_B <= mu_max)
    T_critical = in_firstline_range ? T_of_mu(mu_B) : -Inf

    println("mu_B = $mu_B MeV, in_firstline_range = $in_firstline_range, T_crit = $T_critical")

    X00 = [-1.8, -1.8, -2.2, 0.01, 0.01]
    X00 = Tmu(0.1/hc, 0.0/hc, X00, q, ints)
    phi = X00[1:3]; Phi1 = X00[4]; Phi2 = X00[5]
    P0 = -Omega(phi, Phi1, Phi2, 0.1/hc, 0.0/hc, q, ints)
    println("Reference pressure P0 = $P0 fm^4")

    T_start = Ts[1]
    initial_confined = in_firstline_range ? (T_start < T_critical) : (mu_B < mu_min)
    X0 = initial_confined ? copy(X_CONF) : copy(X_DECONF)
    println("Initial phase: $(initial_confined ? "confined" : "deconfined")")

    data = zeros(length(Ts), 9)  # T, mu_B, Fluctuations...
    crossed = false

    for (i, T) in enumerate(Ts)
        println("T = $T")

        # 检测穿越相变线
        if in_firstline_range && (!crossed) && (T >= T_critical)
            crossed = true
            println("  crossing first-order line at T ≈ $T_critical, switching to deconfined phase")
            X0 .= X_DECONF
        end

        # 求解
        NewX = Tmu(T/hc, mu_B/hc, X0, q, ints)
        vals = Fluctuations(1.0001 * NewX, T/hc, mu_B/hc, P0, q, ints)
        
        data[i, :] = [vals...]

        X0 .= NewX
    end
    df = DataFrame(data, [:T, :mu_B, :P, :chi1, :chi2, :chi3, :chi4, :muBT, :P_T4])
    CSV.write("flu/rep_muB=$(round(mu_B,digits=3))_q=$(q)_flu.csv", df)

end



function main_muBT(scale::Float64 = 1.0)
    """
    固定 mu_B / T 比例，扫描温度 T
    """
    Ts = 50.0:2.0:300.0   # 单位:MeV
    mu_Bs = [scale * T for T in Ts]
    ints = get_nodes(128, nodes2=128)
    data = zeros(length(Ts), 9)  # T, mu_B, Fluctuations...
    qs = [1.001, 1.02, 1.05]
    for q in qs 
        X0 = [-1.9, -1.9, -2.2, 0.0038, 0.0038]  # 初始猜测：禁闭相
        X00 = Tmu(0.1/hc, 0.0/hc, X0, q, ints)
        phi = X00[1:3]; Phi1 = X00[4]; Phi2 = X00[5]
        P0 = -Omega(phi, Phi1, Phi2, 0.1/hc, 0.0/hc, q, ints)
        for (i, (T, mu_B)) in enumerate(zip(Ts, mu_Bs))
            println("T = $T, mu_B = $mu_B")
            NewX = Tmu(T/hc, mu_B/hc, X0, q, ints)
            vals = Fluctuations(1.0001 * NewX, T/hc, mu_B/hc, P0, q, ints)
            data[i, :] = [vals...]  # 使用索引 i 赋值到矩阵的第 i 行
            X0 .= NewX
        end
        df = DataFrame(data, [:T, :mu_B, :P, :chi1, :chi2, :chi3, :chi4, :muBT, :P_T4])
        CSV.write("flu/rep_muBT=$(scale)_q=$(q)_flu.csv", df)
    end

end