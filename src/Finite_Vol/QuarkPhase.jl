using Peaks
using NLsolve
using QuadGK
using DataInterpolations
using DataFrames, CSV
using Plots
# using PlotlyJS
using Roots
using FiniteDifferences
using Optim
using Dates



"""
从高温数据中在区间 [T_min, T_max] 上用 *全局*（非局部）最大值确定 T_c：
T_c := argmax_T d(phi_u)/dT,  T ∈ [T_min, T_max].

- 输入 T, phi_u 可以无序，函数内部按 T 升序整理
- 默认仅考虑 T ≥ 110 MeV；可用关键字 T_min, T_max 调整窗口
- 使用 PCHIP(T, phi_u) 作为平滑插值，再用有限差分求导

返回：T_c（若数据不足则返回 NaN）
"""
function find_Tmu(T, phi_u)
    @assert length(T) == length(phi_u) "T, phi_u 长度不一致"

    # 排序与窗口
    p = sortperm(T)
    T = collect(T[p]); phi_u = collect(phi_u[p])
   
    n = length(T)
    if n < 7 || !(first(T) < last(T))
        return NaN
    end
    T_lo, T_hi = first(T), last(T)

    # 插值 φ(T) 与一阶导 g(T)
    φ  = PCHIPInterpolation(phi_u, T)
    fd = FiniteDifferences.central_fdm(5, 1)
    g(t) = fd(φ, t)

    # 仅在开区间内优化，严格避开端点
    ϵ = max(1e-6*(T_hi - T_lo), 1e-6)
    a, b = T_lo + ϵ, T_hi - ϵ
    if !(a < b)
        return NaN
    end

    # --- 多起点粗细结合：先粗扫，再在局部区间用 Brent 精化 ---
    # 粗扫：均匀采样（不含端点）
    N = max(50, min(400, 5n))
    ts = range(a, b; length=N)
    # 避免端点数值波动，去掉两端一点
    ts = collect(ts)[2:end-1]
    gs = map(g, ts)

    # 取若干最佳种子区间（K 个），每个在相邻网格形成的小区间上做 Brent
    K = min(3, length(ts)-1)
    order = sortperm(gs; rev=true)[1:K]              # 按 g 最大的网格点索引
    candidates = Tuple{Float64,Float64}[]
    for k in order
        # 以网格点 k 为中心的最小邻域区间
        L = max(1, k-1); R = min(length(ts), k+1)
        a_loc = max(a, ts[L]); b_loc = min(b, ts[R])
        if a_loc < b_loc
            push!(candidates, (a_loc, b_loc))
        end
    end
    # 去重
    candidates = unique(candidates)

    # 在每个候选区间做 1D Brent（不比较端点，只取极值点）
    best_T, best_val = NaN, -Inf
    for (aa, bb) in candidates
        # 在 (aa, bb) 内部再缩一点，确保不是端点
        δ = 1e-10*(bb - aa)
        aa′, bb′ = aa + δ, bb - δ
        if aa′ >= bb′
            continue
        end
        res = Optim.optimize(t -> -g(t), aa′, bb′, Brent())
        T_star = Optim.minimizer(res)
        # 忽略数值跑到边界的情况
        if T_star <= aa || T_star >= bb
            continue
        end
        val = g(T_star)
        if val > best_val
            best_val = val
            best_T   = T_star
        end
    end

    return best_T
end



function find_Trho(mu, rho)

    @assert length(mu) == length(rho) "mu, rho 长度不一致"

    # ---- 预处理：按 rho 升序、去重 ----
    p = sortperm(rho)
    rho  = collect(rho[p])
    mu  = collect(mu[p])

    rho_min, rho_max = first(rho), last(rho)
    
    inds_max, vals_max = findmaxima(mu)
    inds_min, vals_min = findminima(mu)
    try
        if isempty(inds_max) || isempty(inds_min)
            error("未找到局部极大/极小，无法做 Maxwell")
        end
    catch e
        println("警告: ", e.msg)
        return (
            rho1 = NaN,
            rho2 = NaN,
            rho_minus = NaN,
            rho_plus = NaN,
            mu_star = NaN,
            area_check = NaN,
            note = "无一阶相变或数据不足"
        )
    end
    i_max = inds_max[argmax(vals_max)]
    i_min = inds_min[argmin(vals_min)]
    mu_plus = mu[i_max]
    mu_minus = mu[i_min]
    rho_plus  = rho[i_max]
    rho_minus = rho[i_min]
    mu_mid = (mu[i_max] + mu[i_min]) / 2

    Mu = PCHIPInterpolation(mu, rho)
    xs = range(rho_min, rho_max; length=4101)
    ys = Mu.(xs)
    function outer_roots_at(mu_c::Real)
        roots = Float64[]
        for i in 1:length(xs)-1
            a = ys[i]   - mu_c
            b = ys[i+1] - mu_c
            if a == 0
                push!(roots, xs[i])
            elseif a*b < 0
                t = a / (a - b)                   # 线性比例
                x = xs[i] + t*(xs[i+1]-xs[i])
                push!(roots, x)
            end
        end
        @assert length(roots) >= 2 "当前水平线下交点不足两处；请检查数据/加密网格"
        sort!(roots)
        return first(roots), last(roots)
    end
    
    
    function A(mu_c::Real)
        try
            rho1, rho2 = outer_roots_at(mu_c)
            val, _ = quadgk(r -> Mu(r) - mu_c, rho1, rho2; rtol=1e-8)
            return val
        catch e
            # 如果找不到两个交点，返回一个极大值，nlsolve 会避开这个点
            return 1e10
        end
    end

    NewX = nlsolve((F, x) -> F[1] = A(x[1]), [mu_mid])
    mu_star = NewX.zero[1]
    rho1, rho2 = outer_roots_at(mu_star)
    area_chk, _ = quadgk(r -> Mu(r) - mu_star, rho1, rho2; rtol=1e-10)
    return (
        rho1       = rho1,
        rho2       = rho2,
        rho_minus  = rho_minus,
        rho_plus   = rho_plus,
        mu_minus  = mu_minus,
        mu_plus   = mu_plus,
        mu_star    = mu_star,
        area_check = area_chk,
    )
end

function main_Trho(;R=30.0, el=30.0)
    println("Time:", Dates.now())
    #path = "../../data/FV/T_rho_B_scan_R=$R.dat"
    path = "../../data/FV/T_rho_B_scan_el=$el.dat"
    df = CSV.read(path, DataFrame)
    T = df.T
    mu = df.mu
    rho = df.rho_B

    Ts = unique(T)
    Ts = sort(Ts; rev=false)
    data = zeros(length(Ts), 8) # T, rho1, rho2, rho_minus, rho_plus, mu_star
    valid_rows = 0
    
    for (i, T_i) in enumerate(Ts)
        println("T = $T_i")
        p = findall(T .== T_i)
        result = find_Trho(mu[p], rho[p])
        
        # 检查结果是否包含 note 字段（说明没有找到局部极大/极小值）
        if haskey(pairs(result), :note)
            println("警告: $(result.note) - 跳过温度 T = $T_i")
            continue
        end
        
        valid_rows += 1
        data[valid_rows, :] .= (
            T_i,
            result.rho1,
            result.rho2,
            result.rho_minus,
            result.rho_plus,
            result.mu_star,
            result.mu_minus,
            result.mu_plus,
        )
    end

    # 只保留有效的行
    data = data[1:valid_rows, :]

    df = DataFrame(
        T        = data[:, 1],
        rho1     = data[:, 2],
        rho2     = data[:, 3],
        rho_minus= data[:, 4],
        rho_plus = data[:, 5],
        mu_star  = data[:, 6],
        mu_minus = data[:, 7],
        mu_plus  = data[:, 8],
    )
    #outpath = "../../data/FV/Trho_MaxwellR=$R.dat"

    outpath = "../../data/FV/Trho_Maxwell_el=$el.dat"
    CSV.write(outpath, df)
    println("结果已保存至 $outpath")
end



function main_Tmu(;R=30.0, el=30.0)
    println("Time:", Dates.now())
    path = "../../data/FV/T_mu_B_scan_el=$el.dat"
    df = CSV.read(path, DataFrame)
    T = df.T
    phi_u = df.phi_u
    mu_B = df.mu_B
    mu_B = unique(mu_B)


    data = zeros(length(mu_B), 2) # mu_B, T_c
    for (i, mu_i) in enumerate(mu_B)
        println("mu_B = $mu_i")
        p = findall(df.mu_B .== mu_i)
        T_i = find_Tmu(T[p], phi_u[p])
        println("T_c = $T_i")
        data[i, :] .= (mu_i, T_i)
    end
    df = DataFrame(
        mu_B = data[:, 1],
        T_c  = data[:, 2],
    )
    outpath = "../../data/FV/Tmu_Tc_el=$el.dat"
    CSV.write(outpath, df)
    println("结果已保存至 $outpath")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main_Tmu()
    #main_Trho()
end