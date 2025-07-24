include("constants.jl")
include("tools.jl")
using ForwardDiff



function Omega(phi, Phi1, Phi2, T,mu_B, ints)
    Mu = [mu_B/3, mu_B/3, mu_B/3]
    Mu = reshape(Mu, 3, 1)
    int1 = ints[1]
    int2 = ints[2]
    intsame1 = ints[3]
    intsame2 = ints[4]

    E1 = Energy(phi, int1)
    E2 = Energy(phi, int2)

    chi = chiral(phi)

    U = calc_U(T, Phi1, Phi2)

    E_minus = E2 .- Mu
    E_plus = E2 .+ Mu

    term1 = int_same1' .* (-6 * E1)
    term_log = log.(AA.(E_minus, T, Phi1, Phi2)) .+ log.(AAbar.(E_plus, T, Phi1, Phi2))
    term2 = -2 .* T .* int_same2' .* term_log
    return chi + U + sum(term1) + sum(term2)
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

# 对T求导
function dOmgea_dT(phi, Phi1, Phi2, T, mu_B, ints)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, x, mu_B, ints), T)
end

# 对T2求导
function dOmgea_dT2(phi, Phi1, Phi2, T, mu_B, ints)
    return ForwardDiff.derivative(x -> dOmgea_dT(phi, Phi1, Phi2, x, mu_B, ints), T)
end


# 对mu_B求导
function dOmgea_dmu_B(phi, Phi1, Phi2, T, mu_B, ints)   
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, T, x, ints), mu_B)
end


# 对mu_B2求导
function dOmgea_dmu_B2(phi, Phi1, Phi2, T, mu_B, ints)   
    return ForwardDiff.derivative(x -> dOmgea_dmu_B(phi, Phi1, Phi2, T, x, ints), mu_B)
end
# 对mu_BT求导

function dOmgea_dmu_BT(phi, Phi1, Phi2, T, mu_B, ints)   
    return ForwardDiff.derivative(x -> dOmgea_dmu_B(phi, Phi1, Phi2, x, mu_B, ints), T)
end

# 对T mu_B求导
function dOmgea_Tmu_B(phi, Phi1, Phi2, T, mu_B, ints)   
    return ForwardDiff.derivative(x -> dOmgea_dT(phi, Phi1, Phi2, T, x, ints), mu_B)
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

function Energy(phi, int_data)
    """使用Julia实现的能量函数"""
    p = int_data[1]  # 假设int_data的第一个元素是动量p
    mass = reshape(Mass(phi), 3, 1)
    return sqrt.(p'.^2 .+ mass.^2)
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

function function_Quark_core(x_array, T, mu_B, int_params)
    phi_u, phi_d, phi_s, Phi1, Phi2 = x_array
    phi = [phi_u, phi_d, phi_s]
    dOmega_phi = dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, int_params)
    dOmega_Phi1 = dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, int_params)
    dOmega_Phi2 = dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, int_params)
    return [dOmega_phi[1], dOmega_phi[2], dOmega_phi[3], dOmega_Phi1, dOmega_Phi2]
end