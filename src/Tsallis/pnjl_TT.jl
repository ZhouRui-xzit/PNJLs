
include("constants.jl")

using ForwardDiff # AD


using NaNMath   # nanlog
using FastGaussQuadrature  # Gauss-Legendre 积分
using MeshGrid # meshgrid
using NLsolve # 非线性方程组求解器
using NaNMath

#  gauss 节点 和LL 节点
function gauleg(a, b, n)
    x, w = gausslegendre(n)
    x_mapped = (b - a) / 2 .* x .+ (b + a) / 2
    w_mapped = (b - a) / 2 .* w
    return x_mapped, w_mapped
end

function get_nodes(p_num, n_num)
    p1, w1 = gauleg(0.0, 40.0, p_num)
    ns = range(0, n_num-1)
    wn = ones(length(ns))


    P1, N1 = meshgrid(p1, ns)
    W1, WN = meshgrid(w1, wn)
    Alpha = ones(n_num)
    for i in 2:n_num
        Alpha[i] = 2
    end
    W = W1 .* WN .* Alpha ./ (2 * pi^2)
    ints = [P1, N1, W]
    return ints
end



function AA(x, T, Phi1, Phi2, q)
    """
    夸克分布归一化分母
    
    
    """
    term1 = exp_q(-x / T, q)
    term2 = exp_q(-2.0 * x / T, q)
    term3 = exp_q(-3.0 * x / T, q)

    result = 1.0 + 3.0 * Phi1 * term1 + 3.0 * Phi2 * term2 + term3
    return result
end

function AAbar(x, T, Phi1, Phi2, q)
    """
    反夸克分布归一化分母
    

    """
    term1 = exp_q(-x / T, q)
    term2 = exp_q(-2.0 * x / T, q)
    term3 = exp_q(-3.0 * x / T, q)

    result = 1.0 + 3.0 * Phi2 * term1 + 3.0 * Phi1 * term2 + term3
    return result
end




@inline function exp_q(x, q)
   
    if x>0
        return (1 .+ (q-1).*x) .^(1/(q-1))
    else
        return (1 .+ (1-q).*x) .^(1/(1-q))
    end
end

@inline function log_q(x, q)
  
    if x>0
        return (x .^ (q-1).- 1) ./ (q-1)
    else
        return (x.^(1-q) .- 1) ./(1-q)
    end
end





function dzeta(x, nodes_finite)
    # nodes_finite 是 [0,1] 上的节点和权重
    
    u, w_finite = nodes_finite
    
    # 指数变换: p = -log(1-u)
    p = -log.(1.0 .- u)
    
    # Jacobian 因子
    jacobian = 1.0 ./ (1.0 .- u)
    
    val = @. ( x * sqrt(1.0 + p^2 / x^2) * (2.0*x*atan(p/x) + p*log(x^2 + p^2)) ) /
          ( 2.0 * (-1.0 + exp(2.0*pi*p)) * sqrt(x^2 + p^2) )
    
    term2 = sum( w_finite .* val .* jacobian )
    term1 = x * (-x + 2.0*(-1.0 + x)*log(x)) / 4.0
    
    return term1 + 2.0 * term2
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
    """计算Polyakov-loop势能"""
    Ta = a0 + a1 * (T0/T) + a2 * (T0/T)^2
    Tb = b3 * (T0 / T)^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    
    log_term = NaNMath.log(value)
    
    U = T^4 * (-1/2 * Ta * Phi2 * Phi1 + Tb * log_term)  # 对数有效势

    return U
end




function Omega_mag(M, qf, eB, nodes)
    """计算单个味夸克的磁化能量"""
    x = M^2 / (2 * abs(qf) * eB)
    term1 = -Nc * (abs(qf) * eB)^2 / (2 * π^2) 
    term2 = dzeta(x, nodes) - 0.5 * (x^2 - x) * NaNMath.log(x) + x^2 / 4
    return 1.0 * term1 * term2
end



function Omega(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    # 1. 提前计算常用值

    P, N, W = gauss_nodes
    G = G_eff(eB)
    chi = chiral(phi, G)
    U = calc_U(T, Phi1, Phi2)
    Masses = Mass(phi, G)
    mu = mu_B/3  # 夸克化学势 
    
    # 2. 预分配结果存储空间，减少动态内存分配
    term1 = 0.0
    term_log_sum = 0.0
    term_mag = 0.0
    # 3. 合并循环，避免创建中间数组
    for i in 1:3

        term1 += -3/(8*pi^2) * (Lambda_f*sqrt(Lambda_f^2 + Masses[i]^2)*(2*Lambda_f^2 + Masses[i]^2) 
        - Masses[i]^4 * NaNMath.log((Lambda_f+sqrt(Lambda_f^2+Masses[i]^2))/Masses[i]))

        # 有限温度
        p_i = @. sqrt(2 * N * qf[i] * eB + P^2)
        E_i = @. sqrt(p_i^2 + Masses[i]^2)

        E_minus = @. E_i - mu
        E_plus = @. E_i + mu
        

        log_sum = @. (log_q(AA(E_minus, T, Phi1, Phi2, q), q) + log_q(AAbar(E_plus, T, Phi1, Phi2, q), q))
        # 把 |q_f| 提到 flavor 循环的外层乘上
        term_log_sum += abs(qf[i]) * sum(log_sum .* W)

        term_mag += Omega_mag(Masses[i], qf[i], eB, zeta_nodes)
    end
    

    term2 = -T * eB * term_log_sum
    #println("chi = $chi, term1 = $term1, term2=$term2, U = $U", ", term_mag = $term_mag")

    return chi + term1 + term2 + U + term_mag
end

# 对phi求导
function dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    return ForwardDiff.gradient(x -> Omega(x, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes), phi)
end


# 对Phi1求导
function dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    return ForwardDiff.derivative(x -> Omega(phi, x, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes), Phi1)
end

# 对Phi2求导
function dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, x, T, mu_B, eB, q, gauss_nodes, zeta_nodes), Phi2)
end


# 对mu_B求导
function dOmgea_dmu_B(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, T, x, eB, q, gauss_nodes, zeta_nodes), mu_B)
end

function function_Quark_core(x_array, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    phi_u, phi_d, phi_s, Phi1, Phi2 = x_array
    phi = [phi_u, phi_d, phi_s]
    dOmega_phi = dOmgea_dphi(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    dOmega_Phi1 = dOmgea_dPhi1(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    dOmega_Phi2 = dOmgea_dPhi2(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    return [dOmega_phi[1], dOmega_phi[2], dOmega_phi[3], dOmega_Phi1, dOmega_Phi2]
end


function Quark_mu(x_array, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    # 动态创建与输入类型兼容的数组
    T_out = promote_type(eltype(x_array), typeof(T), typeof(mu_B))
    fvec = zeros(T_out, 5)

    fvec[1], fvec[2], fvec[3], fvec[4], fvec[5] = function_Quark_core(x_array, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    return fvec
end
 



function Quark_rho(x_array, T, rho_B, eB, q, gauss_nodes, zeta_nodes)
    T_out = promote_type(eltype(x_array), typeof(T), typeof(rho_B))

    fvec = zeros(T_out, 6)
    X0 = x_array[1:5]
    phi = x_array[1:3]
    Phi1 = x_array[4]
    Phi2 = x_array[5]

    mu_B = x_array[6]
    fvec[1], fvec[2], fvec[3], fvec[4], fvec[5] = function_Quark_core(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    rho_now = -1 * dOmgea_dmu_B(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)

    fvec[6] = rho_now / rho0 - rho_B
    return fvec
end

function Tmu(T, mu_B, X0, eB, q, gauss_nodes, zeta_nodes)
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    if res.f_converged == false
        println("Warning: Solver did not converge for T=$(T*hc)")
    end
    NewX = res.zero
    return NewX
end


function Trho(T, rho_B, X0, eB, q, gauss_nodes, zeta_nodes)
    fWrapper(Xs) = Quark_rho(Xs, T, rho_B, eB, q, gauss_nodes, zeta_nodes)
    if rho_B<=1e-8
        res = nlsolve(fWrapper, X0, autodiff=:forward, xtol=1e-10, ftol=1e-13, method=:newton)
        println("xconverged: ", res.x_converged, ", fconverged: ", res.f_converged)
    else 
        res = nlsolve(fWrapper, X0, autodiff=:forward,method=:newton)
    end
    NewX = res.zero
    return NewX
end

