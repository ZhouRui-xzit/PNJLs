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



function main_rep_T(; mu_B::Float64 = 0.0)
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





    for (path, (R, e)) in zip(paths, Res)
        println("Processing firstline file: $path")
        calc_T(mu_B, path, R, e)
    end
end




function calc_T(mu_B::Float64, firstline_path::String, R::Float64, e::Float64)
    """
    固定 mu_B,扫描温度 T
    """
    Ts = 50.0:2.0:300.0
   # Ts2 = 0.1:0.2:1.0
    #Ts = unique(vcat(collect(Ts2), collect(Ts)))
    a, b, c = parametrize_deformation(R, e; para=3.0, scale=-1.0)
    ints = get_nodes_el_hard(128, a, b, c, modes="D")

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

    # 输运系数
    eta_sb = Float64[]; sigma_Tb = Float64[];
    eta_sa = Float64[]; sigma_Ta = Float64[];
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
        vals = Ther_Rep(NewX * 1.0001, T/hc, mu_B/hc, ints, P0)
        P, E, TA, S, rhoB, Cv, cs2 = vals[3], vals[4], vals[5], vals[6], vals[7], vals[8], vals[9]
        trans = trans_eff(X0, T/hc, mu_B/hc, ints)
        s = DTOmega(X0*1.001, T/hc, mu_B/hc, ints)
        eta_s = trans[3]/s 
        sigma_T = trans[4]
        # 存储数据
        if !crossed
            push!(Tb, T); push!(mub, mu_B); push!(Pb, P); push!(Eb, E); push!(TAb, TA);
            push!(Sb, S); push!(rhob, rhoB); push!(cs2b, cs2)
            push!(CVb, Cv); push!(eta_sb, eta_s); push!(sigma_Tb, sigma_T)
        else
            push!(Ta, T); push!(mua, mu_B); push!(Pa, P); push!(Ea, E); push!(TAa, TA);
            push!(Sa, S); push!(rhoa, rhoB); push!(cs2a, cs2)
            push!(CVa, Cv); push!(eta_sa, eta_s); push!(sigma_Ta, sigma_T)
        end

        X0 .= NewX
    end

    # 写文件
    if crossed
        df_before = DataFrame(T=Tb, mu_B=mub, P=Pb, E=Eb, TA=TAb, S=Sb, rhoB=rhob, cs2=cs2b, CV=CVb, eta_s=eta_sb, sigma_T=sigma_Tb)
        df_after  = DataFrame(T=Ta, mu_B=mua, P=Pa, E=Ea, TA=TAa, S=Sa, rhoB=rhoa, cs2=cs2a, CV=CVa, eta_s=eta_sa, sigma_T=sigma_Ta)
        CSV.write("../../data/FV/cs/rep_muB=$(round(mu_B,digits=3))_R=$(R)_e=$(e)_before.csv", df_before)
        CSV.write("../../data/FV/cs/rep_muB=$(round(mu_B,digits=3))_R=$(R)_e=$(e)_after.csv",  df_after)
        println("Saved two files: before & after crossing.")
    else
        df_all = DataFrame(T=Tb, mu_B=mub, P=Pb, E=Eb, TA=TAb, S=Sb, rhoB=rhob, cs2=cs2b, CV=CVb, eta_s=eta_sb, sigma_T=sigma_Tb)
        CSV.write("../../data/FV/cs/rep_muB=$(round(mu_B,digits=3))_R=$(R)_e=$(e)_all.csv", df_all)
        println("No crossing found, saved single file.")
    end
end