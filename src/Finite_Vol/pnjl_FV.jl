include("constants.jl")


using ForwardDiff # AD
using NaNMath   # nanlog
using FastGaussQuadrature  # Gauss-Legendre 积分
using MeshGrid # meshgrid
using NLsolve # 非线性方程组求解器


#  gauss 节点 
function gauleg(a, b, n)
    x, w = gausslegendre(n)
    x_mapped = (b - a) / 2 .* x .+ (b + a) / 2
    w_mapped = (b - a) / 2 .* w
    return x_mapped, w_mapped
end

function surface(p, m)
    return -1/(8*pi) * (1-2/pi * atan(p/m))
end

function total_mean_curvature_ellipsoid(a,b,c; nθ=200, nφ=400)
    x, wθ = gausslegendre(nθ)                 # x ∈ [-1,1] = cosθ
    φ  = range(0, 2π, length=nφ+1)[1:end-1]
    wφ = 2π/nφ
    acc = 0.0
    for (i, xi) in pairs(x)
        s2 = 1 - xi^2                         # sin^2θ
        row = 0.0
        for ϕ in φ
            row += sqrt(a^2*s2*cos(ϕ)^2 + b^2*s2*sin(ϕ)^2 + c^2*xi^2)
        end
        acc += row * wθ[i] * wφ               # 注意：这里没有 sinθ
    end
    return 2 * acc                             # 不是 4！
end


function curvature(p, m)
    return 1/(12*pi^2) * (1 - (3*p)/(2*m) * (pi/2 - atan(p/m)) )
end

#f(x, mass, R) =  1 + (6*pi^2)/(x * R) * surface(x, mass) + (12*pi^2)/(x * R)^2 * curvature(x, mass)

function f_el(x, mass, a, b ,c) # rho_MRE for elipsoid
    V = 4/3 * pi * a * b * c
    S = 4 * pi * ((a^1.6075 * b^1.6075 + a^1.6075 * c^1.6075 + b^1.6075 * c^1.6075)/3)^(1/1.6075)
    C = total_mean_curvature_ellipsoid(a, b, c)
    term1 = 6 * (x^2 * V / (2*pi^2) + x*S * surface(x, mass) + C * curvature(x, mass))
    return term1 / (6 * V)
end


function f_sph(x, mass, R) # rho_MRE for sphere
    V = 4/3 * pi * R^3
    S = 4 * pi * R^2
    C = 8 * pi * R
    term1 =6 * (x^2 * V / (2*pi^2) + x*S * surface(x, mass) + C * curvature(x, mass))
    return term1/(6 * V)
end


function Find_IR_sph(R, masses)
    IR = zeros(3)
    
    for flavor = 1:3
        mass = masses[flavor]
        
        fWrapper(X) = [f_sph(X[1], mass, R)]  # 保证返回向量
        IR0 = [0.1]
        res = nlsolve(fWrapper, IR0, autodiff=:forward)
        IR[flavor] = res.zero[1]
    end
    return IR
end

function Find_IR_el(a, b, c, masses)
    IR = zeros(3)
    for flavor = 1:3
        mass = masses[flavor]
        fWrapper(X) = [f_el(X[1], mass, a, b, c)]  # 保证返回向量
        IR0 = [0.1]
        res = nlsolve(fWrapper, IR0, autodiff=:forward)
        IR[flavor] = res.zero[1]
    end
    return IR
    
end



function get_nodes_sph(num, R ;modes="m")
    if modes == "D"
        masses = alpha_D
    elseif modes == "N"
        masses = alpha_N
    else
        masses = m0
    end
    IR_u, IR_d, IR_s = Find_IR_sph(R, masses)  # 三种夸克的IR截断
    
    # 定义参数组合：每个夸克两种不同的节点配置
    configs = [
        # 真空积分节点 (从IR到Lambda_f)
        (IR_u, Lambda_f, masses[1], "u_vacuum"),  # u夸克真空
        (IR_d, Lambda_f, masses[2], "d_vacuum"),  # d夸克真空
        (IR_s, Lambda_f, masses[3], "s_vacuum"),  # s夸克真空

        # 有限温度积分节点 (从IR到20.0)
        (IR_u, 20.0, masses[1], "u_thermal"),     # u夸克有限温度
        (IR_d, 20.0, masses[2], "d_thermal"),     # d夸克有限温度
        (IR_s, 20.0, masses[3], "s_thermal")      # s夸克有限温度
    ]
    
    # 使用字典来存储结果，便于按夸克类型和积分类型索引
    int_data = Dict()
    
    # 计算每种配置的节点和权重
    for (ir, upper, mass, key) in configs
        p, w = gauleg(ir, upper, num)
        # 计算有限体积修正
        mre = w .* f_sph.(p, mass, R)
        int_data[key] = [p, mre]
    end
    
    # 也可以返回数组格式，保持与原函数兼容
    # 按u_vacuum, u_thermal, d_vacuum, d_thermal, s_vacuum, s_thermal顺序
    return [
        int_data["u_vacuum"], 
        int_data["u_thermal"],
        int_data["d_vacuum"], 
        int_data["d_thermal"],
        int_data["s_vacuum"], 
        int_data["s_thermal"]
    ]
end


function get_nodes_el(num, a, b, c; modes="m")
    if modes == "D"
        masses = alpha_D
    elseif modes == "N"
        masses = alpha_N
    else
        masses = m0
    end
    IR_u, IR_d, IR_s = Find_IR_el(a, b, c, masses)  # 三种夸克的IR截断

    # 定义参数组合：每个夸克两种不同的节点配置
    configs = [
        # 真空积分节点 (从IR到Lambda_f)
        (IR_u, Lambda_f, masses[1], "u_vacuum"),  # u夸克真空
        (IR_d, Lambda_f, masses[2], "d_vacuum"),  # d夸克真空
        (IR_s, Lambda_f, masses[3], "s_vacuum"),  # s夸克真空

        # 有限温度积分节点 (从IR到20.0)
        (IR_u, 20.0, masses[1], "u_thermal"),     # u夸克有限温度
        (IR_d, 20.0, masses[2], "d_thermal"),     # d夸克有限温度
        (IR_s, 20.0, masses[3], "s_thermal")      # s夸克有限温度
    ]
    
    # 使用字典来存储结果，便于按夸克类型和积分类型索引
    int_data = Dict()
    
    # 计算每种配置的节点和权重
    for (ir, upper, mass, key) in configs
        p, w = gauleg(ir, upper, num)
        # 计算有限体积修正
        mre = w .* f_el.(p, mass, a, b, c)
        int_data[key] = [p, mre]
    end
    
    # 也可以返回数组格式，保持与原函数兼容
    # 按u_vacuum, u_thermal, d_vacuum, d_thermal, s_vacuum, s_thermal顺序
    return [
        int_data["u_vacuum"], 
        int_data["u_thermal"],
        int_data["d_vacuum"], 
        int_data["d_thermal"],
        int_data["s_vacuum"], 
        int_data["s_thermal"]
    ]
end





function Omega(phi, Phi1, Phi2, T, mu_B, ints)
    mu = mu_B / 3  # 化学势
    chi = chiral(phi) # 手征相关项
    U = calc_U(T, Phi1, Phi2) # Polyakov-loop 势能项
    Masses = Mass(phi) # 三种夸克有效质量

    Omega_total = chi + U  # 总的热力学势初始值

    for flavor = 1:3
        # 获取对应夸克的节点
        vacuum_idx = 2*flavor - 1    # 真空节点索引：1,3,5
        thermal_idx = 2*flavor       # 热节点索引：2,4,6
        
        # 提取对应的积分节点
        p_vac, w_vac = ints[vacuum_idx] # 当前味道夸克的真空积分节点
        p_therm, w_therm = ints[thermal_idx] # 当前味道夸克的有限温度积分节点
        
        mass = Masses[flavor] # 当前味道夸克的有效质量
        vacuum_contrib = calculate_vacuum_term(p_vac, w_vac, mass) # 真空贡献
        thermal_contrib = calculate_thermal_term(p_therm, w_therm, mass, T, mu, Phi1, Phi2) # 有限温度贡献
        #println("vac:", vacuum_contrib, " ther:", thermal_contrib)

        Omega_total += -2 * Nc * vacuum_contrib - 2 * thermal_contrib # 有限温度没有Nc因子,i.e. Nc因子在胶子场中(AA,AAbar)
    end
    
    return Omega_total

end


function calculate_vacuum_term(p, w, mass)
    E = sqrt.(p.^2 .+ mass^2)
    integrand = w.* E  # 积分测度包含于w中
    return sum(integrand)
end

function calculate_thermal_term(p, w, mass, T, mu, Phi1, Phi2)
    E = sqrt.(p.^2 .+ mass^2)
    E_minus = E .- mu
    E_plus = E .+ mu

    log_sum = log.(AA(E_minus, T, Phi1, Phi2)) .+ log.(AAbar(E_plus, T, Phi1, Phi2))
    integrand = w .* log_sum  # 积分测度包含于w中
    return T * sum(integrand)
end



# 对phi求导
function dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, ints)
    return ForwardDiff.gradient(x -> Omega(x, Phi1, Phi2, T, mu_B, ints), phi)
end

# 对Phi1求导
function dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, ints)
    return ForwardDiff.derivative(x -> Omega(phi, x, Phi2, T, mu_B, ints), Phi1)
end

# 对Phi2求导
function dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, ints)  
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, x, T, mu_B, ints), Phi2)
end


# 对mu_B求导
function dOmgea_dmu_B(phi, Phi1, Phi2, T, mu_B, ints)   
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, T, x, ints), mu_B)
end




function AA(x, T, Phi1, Phi2)
    """
    夸克分布归一化分母（向量化版本）
    
    支持向量输入，计算每个x对应的分母项
    """
    term1 = exp.(-x ./ T)
    term2 = exp.(-2.0 .* x ./ T)
    term3 = exp.(-3.0 .* x ./ T)

    result = 1.0 .+ 3.0 .* Phi1 .* term1 .+ 3.0 .* Phi2 .* term2 .+ term3
    return result
end

function AAbar(x, T, Phi1, Phi2)
    """
    反夸克分布归一化分母（向量化版本）
    
    支持向量输入，计算每个x对应的分母项
    """
    term1 = exp.(-x ./ T)
    term2 = exp.(-2.0 .* x ./ T)
    term3 = exp.(-3.0 .* x ./ T)

    result = 1.0 .+ 3.0 .* Phi2 .* term1 .+ 3.0 .* Phi1 .* term2 .+ term3
    return result
end


function chiral(phi)
    """计算手征相关量
    phi = [phiu, phid, phis]
    """
    term1 = 2 * G_f * sum(phi.^2) - 4 * K_f * phi[1] * phi[2] * phi[3]
    return term1
end

function Mass(phi)
    """计算三种夸克的有效质量"""
    mass_u = m0[1] - 4 * G_f * phi[1] + 2 * K_f * (phi[2] * phi[3])
    mass_d = m0[2] - 4 * G_f * phi[2] + 2 * K_f * (phi[1] * phi[3])
    mass_s = m0[3] - 4 * G_f * phi[3] + 2 * K_f * (phi[1] * phi[2])
    return [mass_u, mass_d, mass_s]
end



function calc_U(T, Phi1, Phi2)
    """计算极化Polyakov-loop势能"""
    Ta = a0 + a1 * (T0/T) + a2 * (T0/T)^2
    Tb = b3 * (T0 / T)^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    
    log_term = NaNMath.log(value)
    
    U = T^4 * (-1/2 * Ta * Phi2 * Phi1 + Tb * log_term)  # 对数有效势

    return U
end

function function_Quark_core(x_array, T, mu_B, int_params)
    # 计算序参量导数
    phi_u, phi_d, phi_s, Phi1, Phi2 = x_array
    phi = [phi_u, phi_d, phi_s]
    dOmega_phi = dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, int_params)
    dOmega_Phi1 = dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, int_params)
    dOmega_Phi2 = dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, int_params)
    return [dOmega_phi[1], dOmega_phi[2], dOmega_phi[3], dOmega_Phi1, dOmega_Phi2]
end


"""
    计算给定T和mu_B下的序参量
    X0 = [phi_u, phi_d, phi_s, Phi1, Phi2] 初始猜测值
"""

function Quark_mu(x_array, T, mu_B, int_params)
    # 动态创建与输入类型兼容的数组
    T_out = promote_type(eltype(x_array), typeof(T), typeof(mu_B))
    fvec = zeros(T_out, 5)

    fvec[1], fvec[2], fvec[3], fvec[4], fvec[5] = function_Quark_core(x_array, T, mu_B, int_params)
    return fvec
end
 

"""
    计算给定T和rho_B下的序参量
    X0 = [phi_u, phi_d, phi_s, Phi1, Phi2, mu_u, mu_d, mu_s] 初始猜测值
"""

function Quark_rho(x_array, T, rho_B,  gauss_nodes)
    T_out = promote_type(eltype(x_array), typeof(T), typeof(rho_B))

    fvec = zeros(T_out, 8)
    X0 = x_array[1:5]
    phi = x_array[1:3]
    Phi1 = x_array[4]
    Phi2 = x_array[5]

    mu_B = x_array[6] * 3
    fvec[1], fvec[2], fvec[3], fvec[4], fvec[5] = function_Quark_core(X0, T, mu_B, gauss_nodes)
    rho_now = -1 * dOmgea_dmu_B(phi, Phi1, Phi2, T, mu_B,  gauss_nodes)

    fvec[6] = x_array[6] - x_array[7]
    fvec[7] = x_array[7] - x_array[8]
    fvec[8] = rho_now / rho0 - rho_B
    return fvec
end

function Tmu(T, mu_B, X0, ints)
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B,  ints)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    NewX = res.zero

    return NewX
end

function Trho(T, rho_B, X0, ints)
    fWrapper(Xs) = Quark_rho(Xs, T, rho_B, ints)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    NewX = res.zero
    return NewX
end
