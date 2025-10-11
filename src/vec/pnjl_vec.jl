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

function get_nodes(p_num)
    p1, w1 = gauleg(0.0,Lambda_f, p_num)
    p2, w2 = gauleg(0.0,20.0, p_num)
    w1 = w1 .* p1.^2 ./ (2*pi^2)
    w2 = w2 .* p2.^2 ./ (2*pi^2)
    int1 = (p1, w1)
    int2 = (p2, w2)
    return int1, int2
end





function Omega(phi, Phi1, Phi2, rhov, T, mu_B, Gv, ints)
    mu = mu_B / 3  - 2 * Gv * rhov  # 有效化学势
    chi = chiral(phi) # 手征相关项
    U = calc_U(T, Phi1, Phi2) # Polyakov-loop 势能项
    Masses = Mass(phi) # 三种夸克有效质量

    Omega_total = chi + U - Gv * rhov^2  # 总的势能密度
    p1, w1 = ints[1]
    p2, w2 = ints[2]
    for flavor = 1:3
        # 获取对应夸克的节点
        
        mass = Masses[flavor] # 当前味道夸克的有效质量
        vacuum_contrib = calculate_vacuum_term(p1, w1, mass) # 真空贡献
        thermal_contrib = calculate_thermal_term(p2, w2, mass, T, mu, Phi1, Phi2) # 有限温度贡献
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
function dOmgea_dphi(phi, Phi1, Phi2, rhov, T, mu_B, Gv, ints)
    return ForwardDiff.gradient(x -> Omega(x, Phi1, Phi2, rhov, T, mu_B, Gv, ints), phi)
end

# 对Phi1求导
function dOmgea_dPhi1(phi, Phi1, Phi2, rhov, T, mu_B, Gv, ints)
    return ForwardDiff.derivative(x -> Omega(phi, x, Phi2, rhov, T, mu_B, Gv, ints), Phi1)
end

# 对Phi2求导
function dOmgea_dPhi2(phi, Phi1, Phi2, rhov, T, mu_B, Gv, ints)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, x, rhov, T, mu_B, Gv, ints), Phi2)
end


# 对mu_B求导
function dOmgea_dmu_B(phi, Phi1, Phi2, rhov, T, mu_B, Gv, ints)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, rhov, T, x, Gv, ints), mu_B)
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

function function_Quark_core(x_array, rhov, T, mu_B, Gv, int_params)
    # 计算序参量导数
    phi_u, phi_d, phi_s, Phi1, Phi2 = x_array
    phi = [phi_u, phi_d, phi_s]
    dOmega_phi = dOmgea_dphi(phi, Phi1, Phi2, rhov, T, mu_B, Gv, int_params)
    dOmega_Phi1 = dOmgea_dPhi1(phi, Phi1, Phi2, rhov, T, mu_B, Gv, int_params)
    dOmega_Phi2 = dOmgea_dPhi2(phi, Phi1, Phi2, rhov, T, mu_B, Gv, int_params)
    return [dOmega_phi[1], dOmega_phi[2], dOmega_phi[3], dOmega_Phi1, dOmega_Phi2]
end


"""
    计算给定T和mu_B下的序参量
    X0 = [phi_u, phi_d, phi_s, Phi1, Phi2] 初始猜测值
"""

function Quark_mu(x_array, T, mu_B, Gv, int_params)
    # 动态创建与输入类型兼容的数组
    T_out = promote_type(eltype(x_array), typeof(T), typeof(mu_B))
    fvec = zeros(T_out, 6)
    phis = x_array[1:5]
    phi = x_array[1:3]
    Phi1 = x_array[4]
    Phi2 = x_array[5]
    rhov = x_array[6]
    fvec[1:5] = function_Quark_core(phis, rhov, T, mu_B, Gv, int_params)
    rho_now = -1 * dOmgea_dmu_B(phi, Phi1, Phi2, rhov, T, mu_B, Gv, ints)
    fvec[6] = rhov - 3 * rho_now 
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
    rho_now = -1 * dOmgea_dmu_B(phi, Phi1, Phi2, rhov, T, mu_B, Gv, ints)

    fvec[6] = x_array[6] - x_array[7]
    fvec[7] = x_array[7] - x_array[8]
    fvec[8] = rho_now / rho0 - rho_B
    return fvec
end

function Tmu(T, mu_B, X0, Gv, ints)
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B, Gv, ints)
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
