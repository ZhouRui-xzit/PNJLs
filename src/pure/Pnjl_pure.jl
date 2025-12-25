############################################################
# PNJL Omega(T, mu_B) 计算核心（含 minimal 方案）

include("constant.jl")

using FastGaussQuadrature
using ForwardDiff
using NLsolve
using Printf
#  gauss 节点 
function gauleg(a, b, n)
    x, w = gausslegendre(n)
    x_mapped = (b - a) / 2 .* x .+ (b + a) / 2
    w_mapped = (b - a) / 2 .* w
    return x_mapped, w_mapped
end

# ==================== 节点缓存工具 ====================



function get_nodes(p_num; nodes2=100)
    # 真空项积分节点 [0, Lambda_f]
    p1, w1 = gauleg(0.0, Lambda_f, p_num)
    w1 = w1 .* p1.^2 ./ (2*pi^2)

    # 热力学项积分节点 [0, ∞)
    # 使用指数变换: p = -log(1 - x), x ∈ [0, 1]
    x, w_x = gauleg(0.0, 1.0, nodes2)
    
    # 变换后的动量节点
    p2 = -log.(1.0 .- x)
    
    # 完整权重 = w_x × jacobian × p² / (2π²)
    # jacobian = dp/dx = 1/(1-x)
    jacobian = 1.0 ./ (1.0 .- x)
    w2 = w_x .* jacobian .* p2.^2 ./ (2*pi^2)
    
    int1 = (p1, w1)
    int2 = (p2, w2)
    return int1, int2
end



function get_nodes_hard(nodes1, nodes2;IR=20.0)
    # 真空项积分节点 [0, Lambda_f]
    p1, w1 = gauleg(0.0, Lambda_f, nodes1)
    w1 = w1 .* p1.^2 ./ (2*pi^2)

    p2, w2 = gauleg(0.0, IR, nodes2)  # 热力学部分使用 [0,IR] 区间的节点
    w2 = w2 .* p2.^2 ./ (2*pi^2)
       
    int1 = (p1, w1)
    int2 = (p2, w2)
    return int1, int2
end

# =========================== 物理子模块 ===========================
function AA(x, T, Phi1, Phi2)
    # 关键修改：使用 .exp 进行广播
    term1 = exp.(-x ./ T)
    term2 = exp.(-2.0 .* x ./ T)
    term3 = exp.(-3.0 .* x ./ T)
    return 1.0 .+ 3.0 .* Phi1 .* term1 .+ 3.0 .* Phi2 .* term2 .+ term3
end

function AAbar(x, T, Phi1, Phi2)
    # 关键修改：使用 .exp 进行广播
    term1 = exp.(-x ./ T)
    term2 = exp.(-2.0 .* x ./ T)
    term3 = exp.(-3.0 .* x ./ T)
    return 1.0 .+ 3.0 .* Phi2 .* term1 .+ 3.0 .* Phi1 .* term2 .+ term3
end



function chiral(phi)
    # phi = [phi_u, phi_d, phi_s]
    return 2 * G_f * sum(phi.^2) - 4 * K_f * phi[1] * phi[2] * phi[3]
end

function Mass(phi)
    # 三味有效质量（fm^-1）
    mass_u = m0[1] - 4 * G_f * phi[1] + 2 * K_f * (phi[2] * phi[3])
    mass_d = m0[2] - 4 * G_f * phi[2] + 2 * K_f * (phi[1] * phi[3])
    mass_s = m0[3] - 4 * G_f * phi[3] + 2 * K_f * (phi[1] * phi[2])
    return [mass_u, mass_d, mass_s]
end

E_of_p(p, mass) = sqrt(p*p + mass*mass)

function calc_U(T, Phi1, Phi2)
    # 对数型 Polyakov 势
    Ta = a0 + a1 * (T0/T) + a2 * (T0/T)^2
    Tb = b3 * (T0 / T)^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    U = T^4 * (-0.5 * Ta * Phi2 * Phi1 + Tb * log(value))
    return U
end

function calculate_vacuum_term(p, w, mass)
    E = sqrt.(p.^2 .+ mass^2)
    integrand = w.* E  # 积分测度包含于w中
    return sum(integrand)
end



function calculate_thermal_term(p, w, mass, T, mu, Phi1, Phi2)
    # p 已经是变换后的动量节点 [0, ∞)
    # w 已经包含了所有因子: Jacobian × p² / (2π²)
    
    # 计算能量
    E = sqrt.(p.^2 .+ mass^2)
    E_minus = E .- mu
    E_plus = E .+ mu

    # 计算被积函数
    log_sum = log.(AA(E_minus, T, Phi1, Phi2)) .+ log.(AAbar(E_plus, T, Phi1, Phi2))
    integrand = w .* log_sum
    
    return T * sum(integrand)
end


function Omega(phi, Phi1, Phi2, T, mu_B, ints)
    mu = mu_B / 3  # 化学势
    chi = chiral(phi) # 手征相关项
    U = calc_U(T, Phi1, Phi2) # Polyakov-loop 势能项
    Masses = Mass(phi) # 三种夸克有效质量

    Omega_total = chi + U  # 总的热力学势初始值
    p1, w1 = ints[1]  # vac nodes
    p2, w2 = ints[2]  # therm nodes

    for flavor = 1:3

        mass = Masses[flavor] # 当前味道夸克的有效质量
        vacuum_contrib = calculate_vacuum_term(p1, w1, mass) # 真空贡献
        thermal_contrib = calculate_thermal_term(p2, w2, mass, T, mu, Phi1, Phi2)
        
     
        Omega_total += -2 * Nc * vacuum_contrib - 2 * thermal_contrib # 有限温度没有Nc因子,i.e. Nc因子在胶子场中(AA,AAbar)
    end
    
    return Omega_total

end


# =========================== dOmega 导数 ===========================
function dOmega_dphi(phi, Phi1, Phi2, T, mu_B, ints)
    return ForwardDiff.gradient(x -> Omega(x, Phi1, Phi2, T, mu_B, ints), phi)
end

function dOmega_dPhi1(phi, Phi1, Phi2, T, mu_B, ints)
    return ForwardDiff.derivative(x -> Omega(phi, x, Phi2, T, mu_B, ints), Phi1)
end

function dOmega_dPhi2(phi, Phi1, Phi2, T, mu_B, ints)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, x, T, mu_B, ints), Phi2)
end



function dOmega_dmu_B(phi, Phi1, Phi2, T, mu_B, ints)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, T, x, ints), mu_B)
end

# =========================== 方程与求解 ===========================
function function_Quark_core(x_array, T, mu_B, ints)
    # 核心 5 维驻点条件
    phi_u, phi_d, phi_s, Phi1, Phi2 = x_array
    phi = [phi_u, phi_d, phi_s]
    dphi  = dOmega_dphi(phi, Phi1, Phi2, T, mu_B, ints)
    dP1   = dOmega_dPhi1(phi, Phi1, Phi2, T, mu_B, ints)
    dP2   = dOmega_dPhi2(phi, Phi1, Phi2, T, mu_B, ints)
    return [dphi[1], dphi[2], dphi[3], dP1, dP2]
end

function Quark_mu(x_array, T, mu_B, ints)
    # T-mu_B 模式：5 维非线性方程
    T_out = promote_type(eltype(x_array), typeof(T), typeof(mu_B))
    fvec = zeros(T_out, 5)
    fvec[1], fvec[2], fvec[3], fvec[4], fvec[5] = function_Quark_core(x_array, T, mu_B, ints)
    return fvec
end

function Quark_rho(x_array, T, rho_B, ints)
    # T-rho_B 模式：8 维非线性方程（含密度约束）
    T_out = promote_type(eltype(x_array), typeof(T), typeof(rho_B))
    fvec = zeros(T_out, 6)

    X0   = x_array[1:5]
    phi  = x_array[1:3]
    Phi1 = x_array[4]
    Phi2 = x_array[5]

    mu_B = x_array[6]
    fvec[1:5] .= function_Quark_core(X0, T, mu_B, ints)

    rho_now = -1 .* dOmega_dmu_B(phi, Phi1, Phi2, T, mu_B, ints)

    fvec[6] = rho_now / rho0 - rho_B
    return fvec
end

function Tmu(T, mu_B, X0, ints)
    # 给定 T, mu_B（均为 fm^-1）求解 5 维驻点
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B, ints)
    res = nlsolve(fWrapper, X0)
    return res.zero
end

function Trho(T, rho_B, X0, ints)
    # 给定 T, rho_B 求解 8 维驻点
    fWrapper(Xs) = Quark_rho(Xs, T, rho_B, ints)
    res = nlsolve(fWrapper, X0)
    return res.zero
end




