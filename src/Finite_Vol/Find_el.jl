using Peaks
using DataFrames, CSV
using Dates  # 记录时间戳
using Revise

includet("pnjl_FV.jl")

"""
判断在给定温度 T（单位：MeV）、给定 theta 下是否存在 S 型曲线
返回命名元组：(has_s_shape, mu, rho, extrema_count)

注意：
- 外部扫描区间里的 T 使用 MeV；在调用 Trho 前会自动做 T/hc 转换为 fm^-1。
- 函数内部对 rho_range 顺序为从大到小扫描，并在末尾统一按 rho 升序排序。
- mu 这里保存的是 X0[6]*hc（即 μ_B/3 的 MeV 值），保持与你原始脚本一致。
"""
function check_s_shape(T, ints; rho_range=3.00:-0.01:0.01)
    
    X0 = [-1.8, -1.8, -2.2, 0.038, 0.038, 1200/hc]
    mu_vals  = Float64[]
    rho_vals = Float64[]

    for rhoB in rho_range
        try
            # 关键：Trho 需要 fm^-1 的温度；外部传入的 T 为 MeV，这里做 T/hc
            X0 = Trho(T/hc, rhoB, X0, ints)
            # 记录 μ_q ≡ μ_B/3 的 MeV 值（与你原脚本一致）
            push!(mu_vals, X0[6] * hc)
            push!(rho_vals, rhoB)
        catch e
            println("T=$(T) MeV, rho=$(rhoB) 计算失败: $e")
            # 失败时重置初值，避免错误扩散
           X0 = [-1.8, -1.8, -2.2, 0.038, 0.038, 0.038, 1200/hc]
        end
    end

    # 统一按 rho 升序排序，便于极值检测和绘图
    p = sortperm(rho_vals)
    rho_vals = rho_vals[p]
    mu_vals  = mu_vals[p]

    # 检测极大/极小点的个数（S 型存在需同时有极大与极小）
    inds_max, _ = findmaxima(mu_vals)
    inds_min, _ = findminima(mu_vals)
    #println("inds_max: ", mu_vals[inds_max])
    #println("inds_min: ", mu_vals[inds_min])
    extrema_count = length(inds_max) + length(inds_min)
    has_s_shape   = !isempty(inds_max) && !isempty(inds_min)

    return (has_s_shape=has_s_shape, mu=mu_vals, rho=rho_vals, extrema_count=extrema_count)
end

"""
粗略扫描定位 T_CEP 的大致范围（温度单位：MeV）
- T_min, T_max, step 统一以 MeV 给定
- 内部调用 check_s_shape(T; theta=...) 时会自动做 T/hc → fm^-1 的转换
"""
function coarse_scan(T_min, T_max, ints; step=1.0)
    #results = Vector{NamedTuple{(:T,:has_s,:extrema),Tuple{Float64,Bool,Int}}}()
    T_range = T_max:-step:T_min

    println("------------------------------------------------")
    println("开始粗略扫描 T ∈ [$T_min, $T_max] MeV, 步长 = $step MeV")
    println("------------------------------------------------")

    for T in T_range
        check = check_s_shape(T, ints)
        if !check.has_s_shape
            println("T = $(T) MeV, S型: false")
        else check.has_s_shape
            println("T = $(T) MeV, S型: true")
            return (T_lower=T, T_upper=T+1)
        end
    end
    println("未找到存在 S 型曲线的温度区间")
    println("I'm going to lower T_min and try again.")
    T_new_rang = T_min - 1:-step:5.0
    for T in T_new_rang
        check = check_s_shape(T, ints)
        if !check.has_s_shape
            println("T = $(T) MeV, S型: false")
        else check.has_s_shape
            println("T = $(T) MeV, S型: true")
            return (T_lower=T, T_upper=T+1)
        end
    end
    
end

"""
二分精确搜索 T_CEP（温度单位：MeV）
- T_lower < T_upper，为粗略扫描得到的包围区间（单位 MeV）
- tolerance 为最终区间宽度容差（单位 MeV）
- rho_density 为 ρ 的采样密度（无量纲）
- theta 以弧度给出
"""
function binary_search_CEP(T_lower, T_upper,ints; tolerance=0.01, rho_density=0.005)
    @assert T_lower < T_upper "T_lower 必须小于 T_upper"

    iter_results = Vector{NamedTuple{(:T,:has_s,:iter),Tuple{Float64,Bool,Int}}}()
    iter, max_iter = 0, 100

    println("\n------------------------------------------------")
    println("开始二分搜索 T_CEP ∈ [$T_lower, $T_upper] MeV")
    println("目标精度: ±$(tolerance/2) MeV")
    println("-----------------------------------------")

    # 边界验证：下边界必须有 S 型，上边界必须无 S 型
    lower_result = check_s_shape(T_lower, ints; rho_range=3.00:-0.01:0.01)
    upper_result = check_s_shape(T_upper, ints; rho_range=3.00:-0.01:0.01)

    if !lower_result.has_s_shape
        println("错误: 下边界温度 $(T_lower) MeV 不存在 S 型曲线")
        return (T_CEP=NaN, note="下边界温度不存在 S 型曲线")
    end
    if upper_result.has_s_shape
        println("错误: 上边界温度 $(T_upper) MeV 仍存在 S 型曲线")
        return (T_CEP=NaN, note="上边界温度仍存在 S 型曲线")
    end

    println("边界验证通过：T=$(T_lower) MeV (S型=true), T=$(T_upper) MeV (S型=false)")

    while T_upper - T_lower > tolerance && iter < max_iter
        iter += 1
        T_mid = (T_lower + T_upper) / 2

        result = check_s_shape(T_mid, ints; rho_range=3.00:-0.01:0.01)
        push!(iter_results, (T=T_mid, has_s=result.has_s_shape, iter=iter))

        println("迭代 $iter: T = $(T_mid) MeV, S型: $(result.has_s_shape), 区间: [$(T_lower), $(T_upper)] MeV")

        if result.has_s_shape
            # 中点仍有 S 型，说明 CEP 在更高温度一侧
            T_lower = T_mid
        else
            # 中点无 S 型，说明 CEP 在更低温度一侧
            T_upper = T_mid
        end
    end

    println("\n二分搜索完成，共 $iter 次迭代")
    println("最终区间: [$(T_lower), $(T_upper)] MeV, 宽度: $(T_upper - T_lower) MeV")

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
                     paras=paras,
                     coarse_step=1.0,
                     fine_tolerance=0.01,
                     rho_density=0.005)

    # 1) 粗略扫描（T 单位：MeV）
    a, b, c = paras #椭球参数
    ints = get_nodes_el(300, a, b, c;modes="D")
    println("=== 开始粗略扫描，步长: $(coarse_step) MeV, a,b,c=$(paras) ===")
    coarse_result = coarse_scan(T_min, T_max,  ints; step=coarse_step)

    if isnan(coarse_result.T_lower)
        println(coarse_result.note)
        return coarse_result
    end
    println("粗略扫描结果: T_CEP ∈ [$(coarse_result.T_lower), $(coarse_result.T_upper)] MeV @ a,b,c=$(paras)")

    # 2) 二分精确搜索
    println("\n=== 开始精确搜索，目标精度: ±$(fine_tolerance/2) MeV, a,b,c=$(paras) ===")
    fine_result = binary_search_CEP(coarse_result.T_lower, coarse_result.T_upper,ints;
                                    tolerance=fine_tolerance, rho_density=rho_density)

    println("\n最终结果:")
    println("T_CEP = $(fine_result.T_CEP) ± $(fine_result.precision/2) MeV")
    println("迭代次数: $(fine_result.iterations)")

    return fine_result
end

"""
脚本主函数：
- ARGS[1] = 目标精度（单位 MeV，可选，默认 0.01）
- ARGS[2] = theta（弧度，可选，默认 0.0）
- 结果写入 ../data/pure/CEP_result_theta.txt（追加），每次记录 theta 与 T_CEP
"""
function main(R, e;T_min=125.0, T_max=135.0)
    paras = parametrize_deformation(R, e;para=3.0,scale=-1.0)
    tolerance = length(ARGS) > 0 ? parse(Float64, ARGS[1]) : 0.01
  
   
    println("====================================================")
    println(" 临界端点 (CEP) 精确搜索")
    println(" 执行时间: $(Dates.now())")
    println(" 单位约定: 扫描温度使用 MeV；Trho 调用内部自动做 T/hc → fm^-1。")
    println(" 参数: tolerance=$(tolerance) MeV, a,b,c=$(paras)")
    println("====================================================")

    result = find_T_CEP(
        T_min=T_min,    # MeV
        T_max=T_max,    # MeV
        paras=paras,
        coarse_step=1.0,          # MeV
        fine_tolerance=tolerance, # MeV
        rho_density=0.005
    )


    outpath = "../../data/FV/CEP_result_el.csv"
    timestamp = string(Dates.now())  # 或 Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    df_new = DataFrame(R=R, e=e,
                    a=paras[1], b=paras[2], c=paras[3],
                    T_CEP=result.T_CEP, timestamp=timestamp)


    if isfile(outpath)
        df_old = CSV.read(outpath, DataFrame)
        df_all = vcat(df_old, df_new)
        CSV.write(outpath, df_all)
    else
        CSV.write(outpath, df_new)
    end
    println("结果已保存至 $outpath")
    
    return result
end

function single_CEP()
    R = 100.0
    delta = 0.7
    #paras = parametrize_deformation(R, delta;para=3.0,scale=-1.0)
    main(R, delta;T_min=80, T_max=130)
end




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



function unique_CEP()
    """
按 (R, e) 去重 CEP 结果文件 ../../data/FV/CEP_result_el.csv，保留最后一次出现的记录。
    """
    df = CSV.read("../../data/FV/CEP_result_el.csv", DataFrame)

    # 若需按 (R,e) 去重并保留最后一次出现：
    if nrow(df) > 0
        # 反转表，使得 unique 保留原表中“最后出现”的那一行
        df_rev = df[end:-1:1, :]
        df_unique_rev = unique(df_rev, [:R, :e])   # unique 保留第一次（此处为原表的最后一次）
        df = df_unique_rev[end:-1:1, :]
    end
    CSV.write("../../data/FV/CEP_result_el.csv", df)
    println("已按 (R, e) 去重并保留最后一次出现的记录。")
end

function CEP_Rs()
    Rs1 = 100.0:-5.0:35.0
    Rs2 = 35.0:-1.0:15.0
    Rs3 = 15.0:-0.1:14.0
    Rs = unique(vcat(Rs1, Rs2, Rs3))
    delta = 0.0
    T_min=4.0
    T_max=132.0
    for R in Rs
        res = main(R, delta;T_min=T_min, T_max=T_max)
        T_max = res.T_CEP + 1.0
        println("R=$R fm, T_CEP=$(res.T_CEP) MeV") 
    end
    unique_CEP()
end

function CEP_deltas()
    #R = 30.0 # T=110 MeV
    #R = 20 # T= 85.xxx MeV
    R = 20.0
    deltas = 0.0:0.02:1.0
    T_min=4.0
    T_max=90.0
    for delta in deltas
        res = main(R, delta;T_min=T_min, T_max=T_max)
        T_max = res.T_CEP + 1.0
        println("delta=$delta , T_CEP=$(res.T_CEP) MeV")
    end
    unique_CEP()
end