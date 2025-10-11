using Peaks
using DataFrames, CSV
using Dates  # 记录时间戳

include("pnjl_axion.jl")

"""
判断在给定温度 T（单位：MeV）、给定 theta 下是否存在 S 型曲线
返回命名元组：(has_s_shape, mu, rho, extrema_count)

注意：
- 外部扫描区间里的 T 使用 MeV；在调用 Trho 前会自动做 T/hc 转换为 fm^-1。
- 函数内部对 rho_range 顺序为从大到小扫描，并在末尾统一按 rho 升序排序。
- mu 这里保存的是 X0[6]*hc（即 μ_B/3 的 MeV 值），保持与你原始脚本一致。
"""
function check_s_shape(T; theta=0.0, rho_range=3.00:-0.01:0.01)
    # 初始猜测值（内部单位要求：X0[9:11] 为 μ_B/3 的 fm^-1）
    # [sigma_u, sigma_d, sigma_s, eta_u, eta_d, eta_s, Phi1, Phi2, muB/3, aux1, aux2]
    #X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1, 320/hc, 320/hc, 320/hc]
    X0 = [1.61781, 1.61781, 2.22585, -0.87333, -0.87333, -0.03414, 0.00383, 0.00383, 300/hc, 300/hc, 300/hc]
    mu_vals  = Float64[]
    rho_vals = Float64[]

    for rhoB in rho_range
        try
            # 关键：Trho 需要 fm^-1 的温度；外部传入的 T 为 MeV，这里做 T/hc
            X0 = Trho(T/hc, rhoB, X0, theta)
            # 记录 μ_q ≡ μ_B/3 的 MeV 值（与你原脚本一致）
            push!(mu_vals, X0[9] * hc)
            push!(rho_vals, rhoB)
        catch e
            println("T=$(T) MeV, rho=$(rhoB) 计算失败: $e")
            # 失败时重置初值，避免错误扩散
            #X0 = [1.8, 1.8, 2.2, 0.01, 0.01, 0.01, 0.1, 0.1, 320/hc, 320/hc, 320/hc]
            X0 = [1.61781, 1.61781, 2.22585, -0.87333, -0.87333, -0.03414, 0.00383, 0.00383, 300/hc, 300/hc, 300/hc]
        end
    end

    # 统一按 rho 升序排序，便于极值检测和绘图
    p = sortperm(rho_vals)
    rho_vals = rho_vals[p]
    mu_vals  = mu_vals[p]

    # 检测极大/极小点的个数（S 型存在需同时有极大与极小）
    inds_max, _ = findmaxima(mu_vals)
    inds_min, _ = findminima(mu_vals)

    extrema_count = length(inds_max) + length(inds_min)
    has_s_shape   = !isempty(inds_max) && !isempty(inds_min)

    return (has_s_shape=has_s_shape, mu=mu_vals, rho=rho_vals, extrema_count=extrema_count)
end

"""
粗略扫描定位 T_CEP 的大致范围（温度单位：MeV）
- T_min, T_max, step 统一以 MeV 给定
- 内部调用 check_s_shape(T; theta=...) 时会自动做 T/hc → fm^-1 的转换
"""
function coarse_scan(T_min, T_max; theta=0.0, step=1.0)
    results = Vector{NamedTuple{(:T,:has_s,:extrema),Tuple{Float64,Bool,Int}}}()
    T_range = T_max:-step:T_min

    println("------------------------------------------------")
    println("开始粗略扫描 T ∈ [$T_min, $T_max] MeV, 步长 = $step MeV, theta = $theta")
    println("------------------------------------------------")

    for T in T_range
        result = check_s_shape(T; theta=theta)
        push!(results, (T=T, has_s=result.has_s_shape, extrema=result.extrema_count))
        println("T = $(T) MeV, theta = $theta, S型: $(result.has_s_shape), 极值数: $(result.extrema_count)")
    end

    # 寻找 S 型从“有”到“无”的温度区间
    df = DataFrame(results)
    sort!(df, :T)

    for i in 1:length(df.T)-1
        if df.has_s[i] && !df.has_s[i+1]
            println("\n找到 S 型曲线消失区间: [$(df.T[i]) MeV, $(df.T[i+1]) MeV] @ theta=$theta")
            return (T_lower=df.T[i], T_upper=df.T[i+1], df=df)
        end
    end

    println("\n警告: 未找到 S 型曲线消失区间，请调整温度范围 (theta=$theta)")
    return (T_lower=NaN, T_upper=NaN, df=df, note="未找到 S 型曲线消失区间，请调整温度范围")
end

"""
二分精确搜索 T_CEP（温度单位：MeV）
- T_lower < T_upper，为粗略扫描得到的包围区间（单位 MeV）
- tolerance 为最终区间宽度容差（单位 MeV）
- rho_density 为 ρ 的采样密度（无量纲）
- theta 以弧度给出
"""
function binary_search_CEP(T_lower, T_upper; theta=0.0, tolerance=0.01, rho_density=0.005)
    @assert T_lower < T_upper "T_lower 必须小于 T_upper"

    iter_results = Vector{NamedTuple{(:T,:has_s,:iter),Tuple{Float64,Bool,Int}}}()
    iter, max_iter = 0, 100

    println("\n------------------------------------------------")
    println("开始二分搜索 T_CEP ∈ [$T_lower, $T_upper] MeV @ theta=$theta")
    println("目标精度: ±$(tolerance/2) MeV")
    println("------------------------------------------------")

    # 边界验证：下边界必须有 S 型，上边界必须无 S 型
    lower_result = check_s_shape(T_lower; theta=theta, rho_range=3.00:-rho_density:0.01)
    upper_result = check_s_shape(T_upper; theta=theta, rho_range=3.00:-rho_density:0.01)

    if !lower_result.has_s_shape
        println("错误: 下边界温度 $(T_lower) MeV 不存在 S 型曲线 (theta=$theta)")
        return (T_CEP=NaN, note="下边界温度不存在 S 型曲线")
    end
    if upper_result.has_s_shape
        println("错误: 上边界温度 $(T_upper) MeV 仍存在 S 型曲线 (theta=$theta)")
        return (T_CEP=NaN, note="上边界温度仍存在 S 型曲线")
    end

    println("边界验证通过：T=$(T_lower) MeV (S型=true), T=$(T_upper) MeV (S型=false) @ theta=$theta")

    while T_upper - T_lower > tolerance && iter < max_iter
        iter += 1
        T_mid = (T_lower + T_upper) / 2

        result = check_s_shape(T_mid; theta=theta, rho_range=3.00:-rho_density:0.01)
        push!(iter_results, (T=T_mid, has_s=result.has_s_shape, iter=iter))

        println("迭代 $iter: T = $(T_mid) MeV, S型: $(result.has_s_shape), 区间: [$(T_lower), $(T_upper)] MeV, theta=$theta")

        if result.has_s_shape
            # 中点仍有 S 型，说明 CEP 在更高温度一侧
            T_lower = T_mid
        else
            # 中点无 S 型，说明 CEP 在更低温度一侧
            T_upper = T_mid
        end
    end

    println("\n二分搜索完成，共 $iter 次迭代")
    println("最终区间: [$(T_lower), $(T_upper)] MeV, 宽度: $(T_upper - T_lower) MeV, theta=$theta")

    return (T_CEP=T_lower, T_lower=T_lower, T_upper=T_upper,
            precision=T_upper - T_lower, iterations=iter,
            df=DataFrame(iter_results))
end

"""
主入口：寻找临界端点（温度单位：MeV）
- T_min, T_max, coarse_step 均以 MeV 指定
- fine_tolerance 为最终 T 区间宽度的目标精度（MeV）
- rho_density 为 ρ 采样步长
- theta 以弧度给出
"""
function find_T_CEP(; T_min=120.0, T_max=140.0,
                     theta=0.0,
                     coarse_step=1.0,
                     fine_tolerance=0.01,
                     rho_density=0.005)

    # 1) 粗略扫描（T 单位：MeV）
    println("=== 开始粗略扫描，步长: $(coarse_step) MeV, theta=$theta ===")
    coarse_result = coarse_scan(T_min, T_max; theta=theta, step=coarse_step)

    if isnan(coarse_result.T_lower)
        println(coarse_result.note)
        return coarse_result
    end
    println("粗略扫描结果: T_CEP ∈ [$(coarse_result.T_lower), $(coarse_result.T_upper)] MeV @ theta=$theta")

    # 2) 二分精确搜索
    println("\n=== 开始精确搜索，目标精度: ±$(fine_tolerance/2) MeV, theta=$theta ===")
    fine_result = binary_search_CEP(coarse_result.T_lower, coarse_result.T_upper;
                                    theta=theta, tolerance=fine_tolerance, rho_density=rho_density)

    println("\n最终结果:")
    println("T_CEP(theta=$(theta)) = $(fine_result.T_CEP) ± $(fine_result.precision/2) MeV")
    println("迭代次数: $(fine_result.iterations)")

    return fine_result
end

"""
脚本主函数：
- ARGS[1] = 目标精度（单位 MeV，可选，默认 0.01）
- ARGS[2] = theta（弧度，可选，默认 0.0）
- 结果写入 ../data/pure/CEP_result_theta.txt（追加），每次记录 theta 与 T_CEP
"""
function main()
    tolerance = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 0.01
    theta     = length(ARGS) > 1 ? parse(Float64, ARGS[2]) : pi/3

    println("====================================================")
    println(" 临界端点 (CEP) 精确搜索")
    println(" 执行时间: $(Dates.now())")
    println(" 单位约定: 扫描温度使用 MeV；Trho 调用内部自动做 T/hc → fm^-1。")
    println(" 参数: tolerance=$(tolerance) MeV, theta=$(theta) (rad)")
    println("====================================================")

    result = find_T_CEP(
        T_min=131.0,    # MeV
        T_max=140.0,    # MeV
        theta=theta,
        coarse_step=1.0,          # MeV
        fine_tolerance=tolerance, # MeV
        rho_density=0.005
    )

    # 目录：../data/pure/
    outdir  = joinpath(@__DIR__, "..", "data", "pure")
    mkpath(outdir)
    outfile = joinpath(outdir, "CEP_result_theta.txt")

    open(outfile, "a") do io
        println(io, "==== CEP Search @ $(Dates.now()) ====")
        println(io, "theta = $(theta) (rad)")
        println(io, "T_CEP(theta) = $(result.T_CEP) ± $(result.precision/2) MeV")
        println(io, "区间: [$(result.T_lower), $(result.T_upper)] MeV")
        println(io, "迭代次数: $(result.iterations)")
        println(io)  # 空行分隔
    end

    println("\n结果已追加保存至 $(outfile)")
end

# 作为脚本执行
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
