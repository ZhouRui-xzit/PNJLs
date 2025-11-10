using Revise
using Dierckx
using DifferentialEquations
using Plots, LaTeXStrings
using CSV, DataFrames

const hc = 197.33
const to_geo = 2.6115e-4 / hc # MeV/fm³ to km⁻² in geometric units
const mass_sun = 1.47664  # Solar mass in km (geometric units)


function logspace(start, stop, n)
    return 10 .^ range(start, stop, length=n)
end

# 初始化 EOS
function init_EOS(en_arr, p_arr)
    # P → E 的三次样条插值
    sort_ind = sortperm(p_arr)
    p_sorted = p_arr[sort_ind]
    en_sorted = en_arr[sort_ind]
    en_dens = Spline1D(p_sorted, en_sorted, k=3)
    
    # E → P 的三次样条插值
    sort_ind2 = sortperm(en_arr)
    en_sorted2 = en_arr[sort_ind2]
    p_sorted2 = p_arr[sort_ind2]
    press = Spline1D(en_sorted2, p_sorted2, k=3)
    
    return en_dens, press
end

# TOV 方程
function tov_eq!(du, u, p, r)
    P, m = u
    en_dens, press = p  # 从参数中获取 B_eff
    
    # 检查压强范围
    if P <= 0 || r <= 0
        du[1] = 0.0
        du[2] = 0.0
        return
    end
    
    # 能量密度
    eden = en_dens(P)
    
    # TOV 方程
    dPdr = -(eden + P) * (m + 4π * r^3 * P) / (r * (r - 2 * m))
    dmdr = 4π * r^2 * eden

    du[1] = dPdr
    du[2] = dmdr
end
#4D-EGB TOV 方程
function EGB_tov_eq!(du, u, p, r)
    P, m = u
    en_dens, press, alpha = p
    # alpha 是 EGB 修正参数
    # 检查压强范围
    if P <= 0 || r <= 0
        du[1] = 0.0
        du[2] = 0.0
        return
    end
    
    # 能量密度
    eden = en_dens(P)
    Gamma = sqrt(1 + 8 * alpha * m / r^3)
    # TOV 方程
    term1 = (eden + P) * (r^3 * (Gamma + 8 * pi * alpha * P - 1) - 2 * alpha * m)
    term2 = r^2 * Gamma * (r^2 * (Gamma - 1) - 2 * alpha)
    dPdr =  term1 / term2
    dmdr = 4 * pi * r^2 * eden

    du[1] = dPdr
    du[2] = dmdr
end



function solve_TOV(c_dens, en_dens, press; rmax=20, rtol=1e-5, dr=0.05, alpha=6.0)
    P0 = press(c_dens)
    eden0 = en_dens(P0)
    m0 = 4π * dr^3 * eden0 / 3.0 
    r_span = (dr, rmax)
    function condition(u, t, integrator)
        return  u[1] - 1e-12  # 夸克星边界条件
    end
    
    function affect!(integrator)
        terminate!(integrator)
    end
    
    cb = ContinuousCallback(condition, affect!)
   
    #prob = ODEProblem(tov_eq!, [P0, m0], r_span, (en_dens, press)) # for Einstein TOV
    prob = ODEProblem(EGB_tov_eq!, [P0, m0], r_span, (en_dens, press, alpha))

    sol = solve(prob, Tsit5(),
               reltol=rtol,
               abstol=1e-10,
               callback=cb, 
               dtmax=dr)
    
    # 直接使用 callback 停止点（不需要事后检查）
    R = sol.t[end] 
    M = sol[2, end] / mass_sun  # 转换为太阳质量单位
    
    return R, M, sol
end

function main(;B_eff=10,alpha=1e-6)
    #B_eff = 70  # 有效袋常数 MeV/fm³ for energy density shift
    #alpha = 1 # EGB 修正参数
    # some alpha : 0.0 1-6   
    # 加载 EOS 数据（单位：MeV/fm³）
    #theta = 0.0
    df = CSV.read("EOS_Rv=0.5_theta=3.141592653.csv", DataFrame)

    en_arr = df.E .+ B_eff  
    p_arr = df.P .- B_eff
    
    println("EOS 范围:")
    println("  能量密度: $(minimum(en_arr)) - $(maximum(en_arr)) MEV/fm³")
    println("  压强:     $(minimum(p_arr)) - $(maximum(p_arr)) MEV/fm³")
    out_file = "data/M_R_theta=pi.csv"
   

    # 绘制 EOS（转回 MeV/fm³ 用于显示）
    p1 = plot(en_arr, 
             p_arr, 
             xlabel=L"E \, (\mathrm{MeV/fm}^3)", 
             ylabel=L"P \, (\mathrm{MeV/fm}^3)", 
             label="", lw=2,
             )
    savefig(p1, "fig/EOS.svg")
    en_arr = en_arr .* to_geo
    p_arr = p_arr .* to_geo



    # 初始化 EOS
    en_dens, press = init_EOS(en_arr, p_arr)
    # 生成密度网格（单位：km-2）
    min_density = en_dens(1e-10)  # P=0 对应的能量密度
    max_density = maximum(en_arr)
    
    # 避免边界问题
    min_log = log10(min_density)
    max_log = log10(max_density)
    
    # 1. 高密度区（靠近最大质量）：密集采样 60%
    dens_high = sort(logspace(max_log - 0.6*(max_log - min_log), max_log, 100), rev=true)
    
    # 2. 中密度区：正常采样 30%
    dens_mid = sort(logspace(max_log - 0.9*(max_log - min_log), 
                             max_log - 0.6*(max_log - min_log), 100), rev=true)
    
    # 3. 低密度区（M→0）：稀疏采样 10%
    dens_low = sort(logspace(min_log, max_log - 0.90*(max_log - min_log), 40), rev=true)
    dens_low_low = sort(logspace(min_log, max_log - 0.999*(max_log - min_log), 10), rev=true)
    dens_low_low_low = sort(logspace(min_log, max_log - 0.99999*(max_log - min_log), 10), rev=true)
    # 合并并去重
    dens = unique(sort([dens_high; dens_mid; dens_low; dens_low_low; dens_low_low_low], rev=true))

    # 过滤掉导致负压强的密度点
    dens = filter(ec -> press(ec) > 1e-20, dens)
    
    # 求解 M-R 关系
    Rs = Float64[]
    Ms = Float64[]
    
    println("\n求解 M-R 关系 (α = $alpha)...")
    for (i, ec) in enumerate(dens)
        try
            R, M, sol = solve_TOV(ec, en_dens, press, 
                                 rmax=20, rtol=1e-5, dr=0.05, alpha=alpha)
            if i ==1
                P = sol[1, :] ./(2.6115e-4 / hc) # to mev/fm³
                m = sol[2, :] ./1.47664 # to m_sun
                df = DataFrame(r=sol.t, P=P, m=m)
                CSV.write("data/TOV_solution.csv", df)
            end
            push!(Rs, R)
            push!(Ms, M)

            ec_MeV = ec / to_geo
            P = press(ec) / to_geo
            println("[$i/$(length(dens))] ec = $(round(ec_MeV, sigdigits=6)) MeV/fm³  Pc = $(round(P, digits=6)) MeV/fm³, " *
                   "R = $(round(R, digits=6)) km, M = $(round(M, digits=6)) M_☉")
        catch e
            ec_MeV = ec / to_geo
            @warn "求解失败 at ρc = $(round(ec_MeV, sigdigits=4)) MeV/fm³: $e"
        end
    end
    # 绘制 M-R 图
    p2 = plot(Rs, Ms, 
            xlimit=(0, 20),
            ylimit=(0, 2.5),
             xlabel=L"R", 
             ylabel=L"Mass/($M_\odot$)", 
             lw=2,
             legend=false
             )
    savefig(p2, "fig/M-R.svg")
    df = DataFrame(R=Rs, M=Ms)
    CSV.write(out_file, df)
    # 找到最大质量
    if !isempty(Ms)
        max_M = maximum(Ms)
        max_idx = argmax(Ms)
        max_R = Rs[max_idx]
        println("\n最大质量星:")
        println("  M_max = $(round(max_M, digits=3)) M_☉")
        println("  R_max = $(round(max_R, digits=2)) km")
    end
    println("M-R has been saved to fig/M-R.svg")
end

