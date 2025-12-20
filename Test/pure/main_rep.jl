
includet("../../src/pure/Rep.jl")
using DataFrames
using CSV
using BenchmarkTools
using Revise
using Dierckx


function main_flu(mu_B=0.0)
    Ts = 300.0:-2.0:1.0
    X0 = [-0.01, -0.01, -0.10, 0.80, 0.80]
    ints = get_nodes(200 ;nodes2=200)
    data = zeros(length(Ts), 7)
    NewX = similar(X0)
    println("mu_B = $mu_B MeV")
    for (i, T) in enumerate(Ts)
        println("T = $T")
        NewX = Tmu(T/hc, mu_B/hc, X0, ints)
        data[i, :] = Fluctuations(NewX*1.001, T/hc, mu_B/hc, ints)
        X0 = NewX
    end

    df = DataFrame(data, [:T, :mu_B, :P, :chimu, :chimu2, :chimu3, :chimu4])
    CSV.write("data_mu=$mu_B.csv", df)
    println("Done!")
end






function calc_T(mu_B::Float64 ;firstline_path="../../data/pure/1st/firstline.csv")
    """
    固定 mu_B,扫描温度 T
    """
    Ts = 50.0:2.0:500.0   # 单位:MeV


    
    ints = get_nodes_hard(10, 64;IR=50.0)

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
    X00 = Tmu(0.1/hc, 0.0/hc, X00, ints)
    phi = X00[1:3]; Phi1 = X00[4]; Phi2 = X00[5]
    P0 = -Omega(phi, Phi1, Phi2, 0.1/hc, 0.0/hc, ints)
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
        NewX = Tmu(T/hc, mu_B/hc, X0, ints)
        vals = Ther_Rep(NewX * 1.0001, T/hc, mu_B/hc, P0,ints,)
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
        CSV.write("../../data/FV/cs/rep_muB=$(round(mu_B,digits=3))_before.csv", df_before)
        CSV.write("../../data/FV/cs/rep_muB=$(round(mu_B,digits=3))_after.csv",  df_after)
        println("Saved two files: before & after crossing.")
    else
        df_all = DataFrame(T=Tb, mu_B=mub, P=Pb, E=Eb, TA=TAb, S=Sb, rhoB=rhob, cs2=cs2b, CV=CVb)
        CSV.write("../../data/FV/cs/rep_muB=$(round(mu_B,digits=3))_all.csv", df_all)
        println("No crossing found, saved single file.")
    end
end