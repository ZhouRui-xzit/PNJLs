using Peaks
using NLsolve
using QuadGK
using DataInterpolations
using DataFrames, CSV
using Roots

"""
识别 S 形区域，只负责：
- 全局按 rho 排序
- 找出 mu(rho) 的所有极大值、极小值
- 按 rho 顺序将第 i 个极大和第 i 个极小配成第 i 个 S
- 返回每个 S 的极大值、极小值信息和全局排序后的数据
不设任何 rho 窗口边界。
"""
function identify_s_shape_regions(rho::AbstractVector, mu::AbstractVector)
    @assert length(rho) == length(mu) "rho, mu 长度不一致"

    # 全局排序
    p = sortperm(rho)
    rho_sorted = rho[p]
    mu_sorted  = mu[p]

    # 全局极值
    inds_max, vals_max = findmaxima(mu_sorted)
    inds_min, vals_min = findminima(mu_sorted)

    n_max = length(inds_max)
    n_min = length(inds_min)

    if n_max < 1 || n_min < 1
        println("  未找到足够的极值点")
        return NamedTuple[]
    end

    n_s = min(n_max, n_min)

    println("  发现 $n_max 个极大值, $n_min 个极小值")
    println("  按 rho 顺序配对得到 $n_s 个 S 形")

    # 按 rho 顺序对极大、极小排序
    max_sorted = sort(collect(zip(inds_max, vals_max)), by = x -> x[1])
    min_sorted = sort(collect(zip(inds_min, vals_min)), by = x -> x[1])

    regions = NamedTuple[]

    for i in 1:n_s
        idx_max, val_max = max_sorted[i]
        idx_min, val_min = min_sorted[i]

        # S 形中心，用于后面 Maxwell 时选根
        rho_center = 0.5 * (rho_sorted[idx_max] + rho_sorted[idx_min])

        println("  S 形 #$i:")
        println("    极大值: 索引 = $idx_max, rho = $(round(rho_sorted[idx_max], digits=3)), mu = $(round(val_max, digits=3))")
        println("    极小值: 索引 = $idx_min, rho = $(round(rho_sorted[idx_min], digits=3)), mu = $(round(val_min, digits=3))")

        region = (
            s_id = i,
            idx_max = idx_max,
            idx_min = idx_min,
            mu_max = val_max,
            mu_min = val_min,
            rho_center = rho_center,
            rho_data = rho_sorted,
            mu_data  = mu_sorted,
            # 只是方便输出查看，不作为边界使用
            rho_range = (
                min(rho_sorted[idx_max], rho_sorted[idx_min]),
                max(rho_sorted[idx_max], rho_sorted[idx_min])
            )
        )

        push!(regions, region)
    end

    return regions
end

function maxwell_single_s(rho_global::AbstractVector,
                          mu_global::AbstractVector,
                          region)

    # 全局排序后的数据
    rho_sorted = region.rho_data
    mu_sorted  = region.mu_data

    idx_max = region.idx_max
    idx_min = region.idx_min
    mu_plus = region.mu_max
    mu_minus = region.mu_min

    # 确保 mu_plus > mu_minus
    if mu_plus < mu_minus
        mu_plus, mu_minus = mu_minus, mu_plus
        idx_max, idx_min = idx_min, idx_max
    end

    # 这一个 S 的极值所在的 rho 位置
    rho_max = rho_sorted[idx_max]
    rho_min = rho_sorted[idx_min]
    rho_left_ext  = min(rho_max, rho_min)
    rho_right_ext = max(rho_max, rho_min)

    amplitude = mu_plus - mu_minus  # 只用来确定搜索区间，不再作为过滤条件

    # 全局 PCHIP 插值 μ(ρ)
    Mu = PCHIPInterpolation(mu_sorted, rho_sorted)
    rho_min_global = first(rho_sorted)
    rho_max_global = last(rho_sorted)

    # 密集网格，用于找交点
    xs = range(rho_min_global, rho_max_global; length = 4001)
    ys = Mu.(xs)

    # 给定 μ_c，在全局上找所有交点，并在极大/极小之外选出左右根
    function outer_roots_at(mu_c::Real)
        roots = Float64[]
        for i in 1:length(xs)-1
            a = ys[i] - mu_c
            b = ys[i+1] - mu_c
            if abs(a) < 1e-12
                push!(roots, xs[i])
            elseif a * b < 0
                t = a / (a - b)
                x = xs[i] + t * (xs[i+1] - xs[i])
                push!(roots, x)
            end
        end

        unique!(roots)
        sort!(roots)

        # 按你说的物理要求：
        # ρ1 必须在 ρ_left_ext 左侧，ρ2 必须在 ρ_right_ext 右侧
        left_candidates  = [x for x in roots if x < rho_left_ext]
        right_candidates = [x for x in roots if x > rho_right_ext]

        if isempty(left_candidates) || isempty(right_candidates)
            error("在 μ=$(round(mu_c, digits=6)) 处无法在极值之外找到左右交点")
        end

        rho1 = maximum(left_candidates)
        rho2 = minimum(right_candidates)

        return rho1, rho2
    end

    # Maxwell 面积函数 A(μ_c) = ∫_{ρ1}^{ρ2} [μ(ρ) - μ_c] dρ
    function A(mu_c::Real)
        rho1, rho2 = outer_roots_at(mu_c)
        val, _ = quadgk(r -> Mu(r) - mu_c, rho1, rho2; rtol = 1e-8, atol = 1e-10)
        return val
    end

    # 在 (μ_min, μ_max) 内搜索 μ*
    # 留一点边界，避免 μ* 贴到极值点本身（理论上 μ*∈(μ_min, μ_max)）
    mu_left  = mu_minus + 0.05 * amplitude
    mu_right = mu_plus  - 0.05 * amplitude

    A_left  = try A(mu_left)  catch; NaN end
    A_right = try A(mu_right) catch; NaN end

    if !isfinite(A_left) || !isfinite(A_right)
        return (
            rho1 = NaN, rho2 = NaN,
            rho_minus = rho_min, rho_plus = rho_max,
            mu_minus = mu_minus, mu_plus = mu_plus,
            mu_star = NaN, area_check = NaN,
            converged = false,
            note = "Maxwell 面积在搜索区间端点不可用"
        )
    end

    # 要求 A(μ) 在 (μ_min, μ_max) 内确实有符号变化
    if A_left * A_right > 0
        return (
            rho1 = NaN, rho2 = NaN,
            rho_minus = rho_min, rho_plus = rho_max,
            mu_minus = mu_minus, mu_plus = mu_plus,
            mu_star = NaN, area_check = NaN,
            converged = false,
            note = "A(μ) 在 (μ_min, μ_max) 内无符号变化根"
        )
    end

    # 用二分法求 μ*，保证 μ* ∈ (μ_min, μ_max)
    mu_star = try
        find_zero(A, (mu_left, mu_right), Bisection(); xtol = 1e-8)
    catch e
        return (
            rho1 = NaN, rho2 = NaN,
            rho_minus = rho_min, rho_plus = rho_max,
            mu_minus = mu_minus, mu_plus = mu_plus,
            mu_star = NaN, area_check = NaN,
            converged = false,
            note = "Maxwell 方程求解失败: $e"
        )
    end

    # 最终的 ρ1, ρ2 和面积检查
    rho1, rho2 = outer_roots_at(mu_star)
    area_chk, _ = quadgk(r -> Mu(r) - mu_star, rho1, rho2; rtol = 1e-10)

    return (
        rho1 = rho1, rho2 = rho2,
        rho_minus = rho_min, rho_plus = rho_max,
        mu_minus = mu_minus, mu_plus = mu_plus,
        mu_star = mu_star, area_check = area_chk,
        converged = true,
        note = ""
    )
end


"""
对包含多个 S 形的 mu-rho 曲线进行 Maxwell 构造。
逻辑：
1. identify_s_shape_regions 只判断 S 的个数和每个 S 的极大/极小；
2. 对每个 S，直接调用 maxwell_single_s；
3. 不在 S 之间做任何 Maxwell 配对，不设 rho 窗口。
"""
function maxwell_multiple_s(rho::AbstractVector, mu::AbstractVector)
    @assert length(rho) == length(mu) "rho, mu 长度不一致"

    regions = identify_s_shape_regions(rho, mu)

    if isempty(regions)
        println("未检测到 S 形区域")
        return NamedTuple[]
    end

    results = NamedTuple[]

    for region in regions
        s_id = region.s_id

        println()
        println("  S 形 #$s_id: 极值区间 rho ≈ [$(round(region.rho_range[1], digits=3)), $(round(region.rho_range[2], digits=3))]")
        println("    mu_max = $(round(region.mu_max, digits=3)), mu_min = $(round(region.mu_min, digits=3))")

        res = maxwell_single_s(rho, mu, region)

        res2 = merge(
            (s_id = s_id,
             rho_range = region.rho_range,
             s_type = length(regions) == 1 ? "single" : "multiple"),
            res
        )

        if res2.converged && isfinite(res2.mu_star)
            println("    Maxwell 成功: mu* = $(round(res2.mu_star, digits=4))")
            println("      共存相 rho ∈ [$(round(res2.rho1, digits=4)), $(round(res2.rho2, digits=4))]")
            println("      面积检查 A(mu*) = $(res2.area_check)")
        else
            println("    Maxwell 失败: $(res2.note)")
        end

        push!(results, res2)
    end

    return results
end



"""
批量处理 T-rho 扫描数据，支持多 S 形 Maxwell 构造
"""
function process_Trho_multiple_maxwell(datafile::String; outfile::String = "")
    df = CSV.read(datafile, DataFrame)

    Ts = unique(df.T)
    sort!(Ts)

    println("="^70)
    println("处理 $(length(Ts)) 个温度点的多 S 形 Maxwell 构造")
    println("="^70)

    all_results = []

    for T_i in Ts
        println("\n" * "="^70)
        println("温度 T = $T_i MeV")
        println("="^70)

        p = findall(df.T .== T_i)
        rho = df.rho_B[p]
        mu  = df.mu_B[p]

        maxwell_results = maxwell_multiple_s(rho, mu)

        for res in maxwell_results
            push!(all_results, merge((T = T_i,), res))
        end
    end

    if isempty(all_results)
        println("\n未找到任何有效的 Maxwell 构造")
        return nothing
    end

    df_out = DataFrame(all_results)

    if isempty(outfile)
        dir = dirname(datafile)
        base = splitext(basename(datafile))[1]
        outfile = joinpath(dir, "$(base)_Maxwell_multi.dat")
    end

    CSV.write(outfile, df_out)

    println("\n" * "="^70)
    println("结果已保存: $outfile")
    println("共处理 $(nrow(df_out)) 个 Maxwell 点 (来自 $(length(Ts)) 个温度)")
    println("="^70)

    return df_out
end

function main_multiple_maxwell()
    datafile = "../../data/magnetic/Trho_eB_0.13_T=10.0.dat"
    results = process_Trho_multiple_maxwell(datafile)
    return results
end