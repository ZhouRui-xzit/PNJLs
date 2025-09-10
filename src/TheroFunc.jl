include("constants.jl")
include("tools.jl")
using ForwardDiff
using FiniteDifferences


function Omega2(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    Mu = [mu_B/3, mu_B/3, mu_B/3]
    Mu = reshape(Mu, 3, 1)
    P, N, W = gauss_nodes

    G = G_eff(eB)

    chi = chiral(phi, G)
    
    U = calc_U(T, Phi1, Phi2)

    Masses = reshape(Mass(phi, G), 3, 1)
    p_array = [(2 .* N .* qf[i] .* eB .+ P.^2).^0.5 for i in 1:3]
    E = [sqrt.(p.^2 .+ Masses[i]^2) for (i, p) in enumerate(p_array)]


    # 计算每种夸克的贡献
    term1_contributions = [qf[i] .* E[i] .* Lambda_f^20 ./ (p_array[i].^20 .+ Lambda_f^20) for i in 1:3]

    # 对每种夸克分别进行积分并求和
    term1_sum = 0.0
    for i in 1:3
        term1_sum += sum(term1_contributions[i] .* W)
    end
    term1 = -Nc * eB * term1_sum


    E_minus = [E[i] .- Mu[i] for i in 1:3]
    E_plus = [E[i] .+ Mu[i] for i in 1:3]

    term_log_contributions = [qf[i] .* (log.(AA(E_minus[i], T, Phi1, Phi2)) .+ log.(AAbar(E_plus[i], T, Phi1, Phi2))) for i in 1:3]

    term_log_sum = 0.0
    for i in 1:3
        term_log_sum += sum(term_log_contributions[i] .* W)
    end
    term2 = -T .* eB .* term_log_sum

    return chi + term1 + term2 + U
end



function Omega(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    # 1. 提前计算常用值
    P, N, W = gauss_nodes
    G = G_eff(eB)
    chi = chiral(phi, G)
    U = calc_U(T, Phi1, Phi2)
    Masses = Mass(phi, G)
    mu = mu_B/3
    
    # 2. 预分配结果存储空间，减少动态内存分配
    term1_sum = 0.0
    term_log_sum = 0.0
    
    # 3. 合并循环，避免创建中间数组
    for i in 1:3
        # 计算动量和能量（一次性计算）
        p_i = @. sqrt(2 * N * qf[i] * eB + P^2)
        E_i = @. sqrt(p_i^2 + Masses[i]^2)
        
        # 4. 使用广播操作代替列表推导
        form_factor = @. Lambda_f^20 / (p_i^20 + Lambda_f^20)
        
        # 5. 直接累加到结果，避免中间数组
        term1_sum += sum(qf[i] .* E_i .* form_factor .* W)
        
        # 计算有限温度贡献
        E_minus = @. E_i - mu
        E_plus = @. E_i + mu
        
        # 使用广播计算对数项
        log_term = @. qf[i] * (log(AA(E_minus, T, Phi1, Phi2)) + 
                               log(AAbar(E_plus, T, Phi1, Phi2)))
        
        term_log_sum += sum(log_term .* W)
    end
    
    # 6. 计算最终结果
    term1 = -Nc * eB * term1_sum
    term2 = -T * eB * term_log_sum
    #println("chi = $chi, term1 = $term1, term2=$term2, U = $U")
    return chi + term1 + term2 + U
end

# 对phi求导
function dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.gradient(x -> Omega(x, Phi1, Phi2, T, mu_B, eB, gauss_nodes), phi)
end


# 对Phi1求导
function dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.derivative(x -> Omega(phi, x, Phi2, T, mu_B, eB, gauss_nodes), Phi1)
end

# 对Phi2求导
function dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, x, T, mu_B, eB, gauss_nodes), Phi2)
end

# 对T求导
function dOmgea_dT(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, x, mu_B, eB, gauss_nodes), T)
end

# 对T2求导
function dOmgea_dT2(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.derivative(x -> dOmgea_dT(phi, Phi1, Phi2, x, mu_B, eB, gauss_nodes), T)
end


# 对mu_B求导
function dOmgea_dmu_B(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, T, x, eB, gauss_nodes), mu_B)
end


# 对mu_B2求导
function dOmgea_dmu_B2(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.derivative(x -> dOmgea_dmu_B(phi, Phi1, Phi2, T, x, eB, gauss_nodes), mu_B)
end
# 对mu_BT求导

function dOmgea_dmu_BT(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.derivative(x -> dOmgea_dmu_B(phi, Phi1, Phi2, x, mu_B, eB, gauss_nodes), T)
end

# 对T mu_B求导
function dOmgea_Tmu_B(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return ForwardDiff.derivative(x -> dOmgea_dT(phi, Phi1, Phi2, T, x, eB, gauss_nodes), mu_B)
end


# 对phi求导
function num_dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    f(x) = Omega(x, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    method = FiniteDifferences.central_fdm(5, 1)
    return FiniteDifferences.grad(method, f, phi)
end

function num_dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    f(x) = Omega(phi, x, Phi2, T, mu_B, eB, gauss_nodes)
    method = FiniteDifferences.central_fdm(5, 1)
    return method(f, Phi1)  # 直接应用差分方法
end

function num_dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    f(x) = Omega(phi, Phi1, x, T, mu_B, eB, gauss_nodes)
    method = FiniteDifferences.central_fdm(5, 1)
    return method(f, Phi2)  # 直接应用差分方法
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

function G_eff(eB)
    xi = eB / Lambda_QCD^2
    return G_f * (1 + a*xi^2+b*xi^3) / (1+ c*xi^2 + d*xi^4)
end

function chiral(phi, G)
    """计算手征相关量
    phi = [phiu, phid, phis]
    """
    term1 = 2 * G * sum(phi.^2) - 4 * K_f * phi[1] * phi[2] * phi[3]
    return term1
end

function Mass(phi, G)
    """计算三种夸克的有效质量"""
    mass_u = m0[1] - 4 * G * phi[1] + 2 * K_f * (phi[2] * phi[3])
    mass_d = m0[2] - 4 * G * phi[2] + 2 * K_f * (phi[1] * phi[3])
    mass_s = m0[3] - 4 * G * phi[3] + 2 * K_f * (phi[1] * phi[2])
    return [mass_u, mass_d, mass_s]
end




function calc_U(T, Phi1, Phi2)
    """计算极化Polyakov-loop势能"""
    Ta = a0 + a1 * (T0/T) + a2 * (T0/T)^2
    Tb = b3 * (T0 / T)^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    
    log_term = log(value)
    
    U = T^4 * (-1/2 * Ta * Phi2 * Phi1 + Tb * log_term)  # 对数有效势

    return U
end

function function_Quark_core(x_array, T, mu_B, eB, gauss_nodes)
    phi_u, phi_d, phi_s, Phi1, Phi2 = x_array
    phi = [phi_u, phi_d, phi_s]
    dOmega_phi = dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    dOmega_Phi1 = dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    dOmega_Phi2 = dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return [dOmega_phi[1], dOmega_phi[2], dOmega_phi[3], dOmega_Phi1, dOmega_Phi2]
end

function num_function_Quark_core(x_array, T, mu_B, eB, gauss_nodes)
    phi_u, phi_d, phi_s, Phi1, Phi2 = x_array
    phi = [phi_u, phi_d, phi_s]
    dOmega_phi = num_dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)[1]
    dOmega_Phi1 = num_dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    dOmega_Phi2 = num_dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, eB, gauss_nodes)
    return [dOmega_phi[1], dOmega_phi[2], dOmega_phi[3], dOmega_Phi1, dOmega_Phi2]
end