using Peaks
using NLsolve
using QuadGK
using DataInterpolations
using DataFrames, CSV
using Plots
# using PlotlyJS
using Roots
using FiniteDifferences



function find_Tmu(T, phi_u)
        # 检查输入数据一致性
    @assert length(T) == length(phi_u) "T, phi_u 长度不一致"

    # 预处理：按 T 升序排序
    p = sortperm(T)
    T = collect(T[p])
    phi_u = collect(phi_u[p])
    
    # 创建插值函数来平滑数据
    phi_interp = PCHIPInterpolation(phi_u, T)
    methods = FiniteDifferences.central_fdm(5, 1)
    dphi_dT = FiniteDifferences.derivative(methods, phi_interp)
    


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
        mu_star    = mu_star,
        area_check = area_chk,
    )
end

function main1()
    df = CSV.read("../../data/pure/T_rho_B_scan1123.dat", DataFrame)
    T = df.T
    mu = df.mu
    rho = df.rho_B

    Ts = unique(T)
    Ts = sort(Ts; rev=false)
    data = zeros(length(Ts), 6) # T, rho1, rho2, rho_minus, rho_plus, mu_star
    for (i, T_i) in enumerate(Ts)
        println("T = $T_i")
        p = findall(T .== T_i)
        result = find_Trho(mu[p], rho[p])
        println(result)
        data[i, :] .= (
            T_i,
            result.rho1,
            result.rho2,
            result.rho_minus,
            result.rho_plus,
            result.mu_star,
        )
    end

    df = DataFrame(
        T        = data[:, 1],
        rho1     = data[:, 2],
        rho2     = data[:, 3],
        rho_minus= data[:, 4],
        rho_plus = data[:, 5],
        mu_star  = data[:, 6],
    )
    CSV.write("../../data/pure/Trho_Maxwell1123.dat", df)
    println("结果已保存至 ../../data/pure/Trho_Maxwell1123.dat")
end

function main2()
    df = CSV.read("mu_rhoB_T130.9.dat", DataFrame)
    mu = df.mu
    rho = df.rho

    println("T = 130.9")
    result = find_Trho(mu, rho)
    println(result)
end
