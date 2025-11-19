using DataFrames
using CSV
using BenchmarkTools
using Revise

includet("../../src/pure/Rep.jl")



function main1(;mu_B=0.0)
    Ts = 300.0:-2.0:50.0
    

    X0 = [-0.01, -0.01, -0.10, 0.80, 0.80]

    data = zeros(length(Ts), 5)
    NewX = similar(X0)
    println("mu_B = $mu_B MeV")
    for (i, T) in enumerate(Ts)
        println("T = $T")
        NewX = Tmu(T, mu_B, X0)
        data[i, :] = Fluctuations(NewX*1.001, T, mu_B)
        X0 = NewX
    end

    df = DataFrame(data, [:T, :mu_B, :chi21, :chi32, :chi42])
    CSV.write("data_mu=$mu_B.csv", df)
    println("Done!")
end


# ===== 工具：读取一阶线，与最近邻索引 =====
using DelimitedFiles

read_firstline(path::AbstractString) = begin
    raw = readdlm(path)             # 文件两列：mu_1stt  T_1stt（单位：MeV）
    mu_1stt = Float64.(raw[:, 1])
    T_1stt  = Float64.(raw[:, 2])
    mu_1stt, T_1stt
end

nearest_index(T::Float64, Ttab::AbstractVector{<:Real}) = begin
    # 最近邻（与 Fortran 的最近邻选择等价，避免过严相等判断）
    idx = findmin(abs.(T .- Ttab))[2]  # 使用[2]访问索引
    Int(idx)
end



function main2(; mu_B::Float64 = 0.0, firstline_path::AbstractString="../../data/pure/1st.txt")
    # mu_BC = 291*3 MEV
    Ts = 50.0:2.0:300.0   # 单位：MeV
        # ===== 两套典型初值（与 Fortran 一致的“相位外形”）=====
    X_CONF   = [-1.9, -1.9, -2.2, 0.0038, 0.0038]
    X_DECONF = [-0.5, -0.5, -1.8, 0.80,   0.80  ]
    P0 = 20.001   # 初始压力值（与 Fortran 一致）
    # 一阶线表：mu_1stt(T), T_1stt (均为 MeV)
    mu_1stt, T_1stt = read_firstline(firstline_path)
    Tmax_1st = maximum(T_1stt)

    # 初始用禁闭样；进入一阶区再切去禁闭样（与 Fortran 相同）
    X0 = copy(X_CONF)

    data = zeros(length(Ts), 5)
    NewX = similar(X0)
    k_1st = 0   # 进入一阶区计数（Fortran 的 k_1st）
    sol = zeros(length(Ts), 6)
    sol[:, 1] .= Ts
    println("mu_B = $mu_B MeV")
    for (i, T) in enumerate(Ts)
        println("T = $T")

        # -------- 一阶区判据与初值切换（Fortran 逻辑复刻） --------
        # 仅当 T 低于一阶线最高温度时，一阶线才存在
        if T < Tmax_1st
            n  = nearest_index(T, T_1stt)
            μ1 = mu_1stt[n]                      # 该温度对应的一阶线 mu_1st
            if mu_B >= μ1                         # 进入一阶区的高 μ 一侧
                k_1st += 1
                if k_1st == 1
                    # 第一次进入：强制切换到去禁闭样初值（Fortran 中的 -0.5 / 0.8 组合）
                    X0 .= X_DECONF
                else
                    # 已在一阶区内：延拓（Fortran 用 f=7 中心点上一解；这里用上一步 NewX）
                    # 无需改动，当前 X0 已在上一轮末尾被设置为 NewX（见下方）
                end
            end
        end
        # ---------------------------------------------------------

        # 牛顿迭代/求解器
        NewX = Tmu(T/hc, mu_B/hc, X0)
        sol[i, 2:end] = NewX
        # 微扰
        data[i, :] = Ther_Rep(NewX*1.0001, T/hc, mu_B/hc, P0)

        # 下一步的延拓初值（Fortran 的 “k>1 then 用上步中心点解” 的 Julia 等价）
        X0 .= NewX
    end

    df = DataFrame(data, [:T, :mu_B, :c_s, :c_rho, :c_s_rho])
    df_sol = DataFrame(sol, [:T, :phiu, :phid, :phis, :Phi, :Phi_bar])
    CSV.write("rep_mu=$(mu_B).csv", df)
    CSV.write("rep_sol_mu=$(mu_B).csv", df_sol)
    println("Done!")
end


function main3(;T=131.0, firstline_path::AbstractString="../../data/pure/1st.txt")
    muBs = 1.0:1.0:1200.0   # 单位：MeV
    #mu_ceps = range(873.0, 874.0, length=20)  # 临界点附近精细扫描

    #muBs = vcat(muBs, mu_ceps...)  # 合并数组
    muBs = sort(muBs)               # 排序
    # ===== 两套典型初值（与 Fortran 一致的"相位外形"）=====
    X_CONF   = [-1.9, -1.9, -2.2, 0.0038, 0.0038]
    X_DECONF = [-0.5, -0.5, -1.8, 0.80,   0.80  ]

    # 一阶线表：mu_1stt(T), T_1stt (均为 MeV)
    mu_1stt, T_1stt = read_firstline(firstline_path)
    
    # 判断当前温度是否超过临界点温度
    T_CEP = maximum(T_1stt)
    is_above_CEP = T > T_CEP
    
    # 初始用禁闭相初值
    X0 = copy(X_CONF)
    
    data = zeros(length(muBs), 3)
    NewX = similar(X0)
    sol = zeros(length(muBs), 6)
    sol[:, 1] .= muBs
    data[:, 1] .= T
    data[:, 2] .= muBs

    k_1st = 0   # 进入一阶区计数
    
    println("T = $T MeV")
    
    # 只有当温度低于临界点温度时才查找一阶相变点
    mu_1st_at_T = 0.0
    if !is_above_CEP
        # 找到给定温度T对应的一阶相变点的化学势
        n_temp = nearest_index(T, T_1stt)
        mu_1st_at_T = mu_1stt[n_temp]
        println("一阶相变点的化学势: mu_B_1st = $mu_1st_at_T MeV")
    else
        println("温度 T = $T MeV 高于临界点温度 T_CEP = $T_CEP MeV，不存在一阶相变")
    end
    
    for (i, mu_B) in enumerate(muBs)
        
        
        # -------- 一阶区判据与初值切换（仅当温度低于临界点温度时） --------
        if !is_above_CEP && mu_B >= mu_1st_at_T  # 进入一阶区的高μ一侧
            k_1st += 1
            if k_1st == 1
                # 第一次进入：强制切换到去禁闭样初值
                X0 .= X_DECONF
                println("  [进入一阶区，切换到去禁闭相初值]")
            end
        end
        # -----------------------------------------
        println("mu_B = $mu_B MeV")
        # 牛顿迭代/求解器
        NewX = Tmu(T/hc, mu_B/hc, X0)
        sol[i, 2:end] = NewX
        phi = NewX[1:3]
        Phi1 = NewX[4]
        Phi2 = NewX[5]
        data[i, 3] = -Omega(phi, Phi1, Phi2, T/hc, mu_B/hc)
        # 下一步的延拓初值
        X0 .= NewX
    end

    df = DataFrame(data, [:T, :mu_B, :P])
    df_sol = DataFrame(sol, [:mu_B, :phiu, :phid, :phis, :Phi, :Phi_bar])
    CSV.write("../../data/pure/FLU_T=$T.csv", df)
    CSV.write("../../data/pure/sol_T=$T.csv", df_sol)
    println("Done!")
end

