using Peaks
using DataFrames, CSV
using Dates
using Revise

includet("pnjl_magnetic.jl")

"""
增强的S形检测，可识别嵌套S形和多个S形
返回: (has_s_shape, mu, rho, s_shape_type, critical_points)
- s_shape_type: "none", "single", "nested", "double"
- critical_points: Dict包含所有极值点信息
"""
function check_s_shape_enhanced(T, eB, ints, zeta_nodes; rho_range=0.01:0.01:4.0)
    rho_range1 = 4.00:-0.01:0.01
    rho_range2 = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3]
    rho_range = sort(vcat(rho_range2, rho_range1), rev=true)

    X0 = [-0.15025, -0.14023, -2.06331, 0.05312, 0.07734, 1200/hc]
    mu_vals = Float64[]
    rho_vals = Float64[]

    for rhoB in rho_range
        try
            X0 = Trho(T/hc, rhoB, X0, eB*(1000/hc)^2, ints, zeta_nodes)
            #println("rhoB = $rhoB, mu = $(X0[6] * hc)")
            push!(mu_vals, X0[6] * hc)
            push!(rho_vals, rhoB)
        catch e
           X0 = [-0.15025, -0.14023, -2.06331, 0.05312, 0.07734, 1200/hc]
        end
    end

    # 按rho排序
    p = sortperm(rho_vals)
    rho_vals = rho_vals[p]
    mu_vals = mu_vals[p]

    # 寻找极值点
    inds_max, vals_max = findmaxima(mu_vals)
    inds_min, vals_min = findminima(mu_vals)
    println("Found maxima at indices: ", vals_max)
    println("Found minima at indices: ", vals_min)
    n_max = length(inds_max)
    n_min = length(inds_min)
    
    # 判断S形类型
    s_shape_type = "none"
    has_s_shape = false
    
    if n_max >= 1 && n_min >= 1
        has_s_shape = true
        
        if n_max == 1 && n_min == 1
            s_shape_type = "single"
        elseif n_max == 2 && n_min == 2
            # 检查是嵌套还是分离的两个S形
            # 嵌套: 极值点交替出现 (max-min-min-max 或 min-max-max-min)
            # 分离: 极值点分组出现 (max-min...gap...max-min)
            all_extrema = sort(vcat(
                [(idx, "max", mu_vals[idx]) for idx in inds_max],
                [(idx, "min", mu_vals[idx]) for idx in inds_min]
            ), by=x->x[1])
            
            pattern = [x[2] for x in all_extrema]
            
            # 检查间隔: 如果相邻极值点间隔较大,可能是分离的
            gaps = diff([x[1] for x in all_extrema])
            max_gap = maximum(gaps)
            avg_gap = sum(gaps) / length(gaps)
            
            if max_gap > 3 * avg_gap  # 存在异常大的间隔
                s_shape_type = "double"
            else
                s_shape_type = "nested"
            end
        elseif n_max >= 2 || n_min >= 2
            s_shape_type = "multiple"  # 更复杂的情况
        end
    end
    
    critical_points = Dict(
        "maxima" => [(rho_vals[i], mu_vals[i]) for i in inds_max],
        "minima" => [(rho_vals[i], mu_vals[i]) for i in inds_min],
        "n_max" => n_max,
        "n_min" => n_min
    )
    
    return (
        has_s_shape=has_s_shape,
        mu=mu_vals,
        rho=rho_vals,
        s_shape_type=s_shape_type,
        critical_points=critical_points,
        extrema_count=n_max + n_min
    )
end

"""
针对不同S形类型的CEP判据
- single -> nested/double: 第一临界点 (小S形出现)
- nested/double -> single: 中间临界点 (小S形消失,大S形仍存在)  
- single -> none: 传统CEP (S形完全消失)
"""
function analyze_cep_type(lower_result, upper_result)
    lower_type = lower_result.s_shape_type
    upper_type = upper_result.s_shape_type
    
    if lower_type == "none" && upper_type in ["single", "nested", "double"]
        return "standard_cep"  # 标准CEP: 无S形 -> 有S形
    elseif lower_type in ["nested", "double"] && upper_type == "single"
        return "nested_cep"  # 嵌套CEP消失点
    elseif lower_type == "single" && upper_type in ["nested", "double"]
        return "nested_onset"  # 嵌套S形出现点
    elseif lower_type == "single" && upper_type == "none"
        return "standard_cep"  # 标准CEP
    else
        return "unknown_transition"
    end
end

"""
增强的粗扫描，记录S形类型变化
"""
function coarse_scan_enhanced(T_min, T_max, eB, ints, zeta_nodes; step=1.0)
    T_range = T_max:-step:T_min
    transitions = []
    
    println("开始增强粗扫描: T ∈ [$T_min, $T_max] MeV, eB = $eB")
    println("="^60)
    
    prev_result = nothing
    
    for T in T_range
        curr_result = check_s_shape_enhanced(T, eB, ints, zeta_nodes)
        
        println("T = $(T) MeV: $(curr_result.s_shape_type), " *
                "极大值=$(curr_result.critical_points["n_max"]), " *
                "极小值=$(curr_result.critical_points["n_min"])")
        
        # 检测S形类型变化
        if !isnothing(prev_result) && prev_result.s_shape_type != curr_result.s_shape_type
            cep_type = analyze_cep_type(prev_result, curr_result)
            push!(transitions, (
                T_lower=T,
                T_upper=T+step,
                type_lower=curr_result.s_shape_type,
                type_upper=prev_result.s_shape_type,
                cep_type=cep_type
            ))
            println("  ⚠ 检测到相变: $(curr_result.s_shape_type) -> " *
                    "$(prev_result.s_shape_type) [$cep_type]")
        end
        
        prev_result = curr_result
    end
    
    println("="^60)
    println("共发现 $(length(transitions)) 个相变点")
    
    return transitions
end

"""
针对特定相变类型的二分搜索
"""
function binary_search_transition(T_lower, T_upper, eB, ints, zeta_nodes;
                                  tolerance=0.01, target_transition=nothing)
    iter, max_iter = 0, 100
    
    println("\n二分搜索: [$T_lower, $T_upper] MeV")
    if !isnothing(target_transition)
        println("目标相变: $(target_transition)")
    end
    
    while T_upper - T_lower > tolerance && iter < max_iter
        iter += 1
        T_mid = (T_lower + T_upper) / 2
        
        result_mid = check_s_shape_enhanced(T_mid, eB, ints, zeta_nodes)
        result_lower = check_s_shape_enhanced(T_lower, eB, ints, zeta_nodes)
        
        println("  迭代 $iter: T_mid=$(T_mid) MeV, type=$(result_mid.s_shape_type)")
        
        # 根据S形类型调整搜索区间
        if result_mid.s_shape_type == result_lower.s_shape_type
            T_lower = T_mid
        else
            T_upper = T_mid
        end
    end
    
    return (T_cep=T_lower, precision=T_upper-T_lower, iterations=iter)
end

"""
主函数: 寻找所有CEP
"""
function find_all_CEPs(; T_min=80.0, T_max=150.0, eB=0.0,
                        coarse_step=1.0, fine_tolerance=0.01)
    
    # 准备积分节点
    ints = get_nodes(256, 20)
    zeta_nodes = gauleg(0, 1, 256)
    
    println("="^60)
    println("磁场PNJL模型多CEP搜索")
    println("磁场: eB = $eB MeV²")
    println("温度范围: [$T_min, $T_max] MeV")
    println("="^60)
    
    # 1. 粗扫描找所有相变
    transitions = coarse_scan_enhanced(T_min, T_max, eB, ints, zeta_nodes; 
                                      step=coarse_step)
    
    if isempty(transitions)
        println("未发现任何相变点")
        return nothing
    end
    
    # 2. 对每个相变进行精确定位
    results = []
    for (i, trans) in enumerate(transitions)
        println("\n" * "="^60)
        println("精确定位第 $i 个相变点: $(trans.cep_type)")
        println("="^60)
        
        result = binary_search_transition(
            trans.T_lower, trans.T_upper, eB, ints, zeta_nodes;
            tolerance=fine_tolerance,
            target_transition=trans.cep_type
        )
        
        push!(results, merge(trans, result))
        
        println("结果: T_CEP = $(result.T_cep) ± $(result.precision/2) MeV")
    end
    
    # 3. 保存结果
    save_results(results, eB)
    
    return results
end

"""
保存所有CEP结果
"""
function save_results(results, eB)
    outdir = joinpath(@__DIR__, "..", "..", "data", "magnetic")
    mkpath(outdir)
    outfile = joinpath(outdir, "CEP_magnetic_eB_$(eB).txt")
    
    open(outfile, "w") do io
        println(io, "Magnetic PNJL CEP Results @ $(Dates.now())")
        println(io, "eB = $eB MeV²")
        println(io, "="^60)
        
        for (i, res) in enumerate(results)
            println(io, "\nCEP #$i: $(res.cep_type)")
            println(io, "  T_CEP = $(res.T_cep) ± $(res.precision/2) MeV")
            println(io, "  Transition: $(res.type_lower) -> $(res.type_upper)")
            println(io, "  Iterations: $(res.iterations)")
        end
    end
    
    println("\n结果已保存至: $outfile")
end

# 使用示例
function main(;eB=0.13, T_min=130.0, T_max=150.0)
    results = find_all_CEPs(
        T_min=T_min,
        T_max=T_max,
        eB=eB,
        coarse_step=1.0,
        fine_tolerance=0.01
    )
    return results
end