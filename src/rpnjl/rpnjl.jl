include("constants.jl")


using ForwardDiff # AD
using NaNMath   # nanlog
using FastGaussQuadrature  # Gauss-Legendre 积分
using MeshGrid # meshgrid
using NLsolve # 非线性方程组求解器
using NaNMath

#  gauss 节点 
function gauleg(a, b, n)
    x, w = gausslegendre(n)
    x_mapped = (b - a) / 2 .* x .+ (b + a) / 2
    w_mapped = (b - a) / 2 .* w
    return x_mapped, w_mapped
end

function get_nodes(p_num)
    p1, w1 = gauleg(0.0, Lambda, p_num)
    p2, w2 = gauleg(0.0, 20.0, p_num)
    w1 = w1 .* p1.^2 ./ (2*pi^2)
    w2 = w2 .* p2.^2 ./ (2*pi^2)
    int1 = (p1, w1)
    int2 = (p2, w2)
    return int1, int2
end





function Omega(orders, mus, T, ints)

    p1, w1 = ints[1]
    p2, w2 = ints[2]

    phi = orders[1:3]
    Phi1 = orders[4]
    Phi2 = orders[5]
    chi = chiral(phi)
    U = calc_U(T, Phi1, Phi2)
    Masses = Mass(phi)
    Omega_q = 0.0
    for flavor = 1:3
        mu = mus[flavor] 
        mass = Masses[flavor]
        vacuum_term =  calculate_vacuum_term(p1, w1, mass)
        thermal_term = calculate_thermal_term(p2, w2, mass, T, mu, Phi1, Phi2)
        Omega_q += vacuum_term + thermal_term
    end

    Omega_total = chi + U + Omega_q 
    return Omega_total      
end


function calculate_vacuum_term(p, w, mass)
    E = sqrt.(p.^2 .+ mass^2)
    integrand = w.* E  # 积分测度包含于w中
    return -2 * Nc * sum(integrand)
end

function calculate_thermal_term(p, w, mass, T, mu, Phi1, Phi2)
    E = sqrt.(p.^2 .+ mass^2)
    E_minus = E .- mu
    E_plus = E .+ mu

    log_sum = log.(AA(E_minus, T, Phi1, Phi2)) .+ log.(AAbar(E_plus, T, Phi1, Phi2))
    integrand = w .* log_sum  # 积分测度包含于w中
    return -2*T * sum(integrand)
end



# 对序参量求导
function dOmega_dorder(orders, mus, T, ints)
    return ForwardDiff.gradient(x -> Omega(x, mus, T, ints), orders)
end
 


# 对mu求导
function dOmega_dmus(orders, mus, T, ints)
    return ForwardDiff.gradient(x -> Omega(orders, x, T, ints), mus)
end

# 对T求导
function dOmega_dT(orders, mus, T, ints)
    return ForwardDiff.derivative(x -> Omega(orders, mus, x, ints), T)
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
    
    term1 = g_S * sum(phi[1:3].^2)
    term2 = -g_D/2 * (phi[1] * phi[2] * phi[3])
    term3 = 3*g1/2 * (sum(phi[1:3].^2))^2
    term4 = 3*g2 * sum(phi[1:3].^4)

    return term1 + term2 + term3 + term4
end




function Mass(phi)
    """计算三种夸克的有效质量"""
    mass_u = mass_u0 - 2*g_S*phi[1]+g_D/4*phi[2]*phi[3]-2*g1*phi[1]*(sum(phi[1:3].^2))-4*g2*phi[1]^3
    mass_d = mass_d0 - 2*g_S*phi[2]+g_D/4*phi[1]*phi[3]-2*g1*phi[2]*(sum(phi[1:3].^2))-4*g2*phi[2]^3
    mass_s = mass_s0 - 2*g_S*phi[3]+g_D/4*phi[1]*phi[2]-2*g1*phi[3]*(sum(phi[1:3].^2))-4*g2*phi[3]^3

    return [mass_u, mass_d, mass_s]
end


function calc_U(T, Phi1, Phi2)
    b2 = a0 + a1 * T0/T * exp(-a2 * T/T0)
    term1 = -b2/2 * Phi1 * Phi2 - b3/6 * (Phi1^3 + Phi2^3) + b4/4 * (Phi1 * Phi2)^2
    J = (27/(24*pi^2)) * (1 - 6*(Phi1*Phi2) + 4*(Phi1^3 + Phi2^3) - 3*(Phi1*Phi2)^2)
    return T^4*term1 - T^4 * kappa * NaNMath.log(J)
end







"""
    计算给定T和mu_B下的序参量
    X0 = [phi_u, phi_d, phi_s, Phi1, Phi2] 初始猜测值
"""

function Quark_mu(X0, mu_B, T, ints)

    T_out = promote_type(eltype(X0), typeof(T), typeof(mu_B))
    orders = X0[1:5]
    mus = [1/3*mu_B, 1/3*mu_B, 1/3*mu_B]

    fvec = zeros(T_out, 5)
    fvec[1:5] = dOmega_dorder(orders, mus, T, ints)

    return fvec
end




"""
    计算给定T和rho_B下的序参量
    X0 = [phi_u, phi_d, phi_s, Phi1, Phi2, muB, muQ] 初始猜测值
"""

function Quark_rho(X0, T, rho, ints)
    T_out = promote_type(eltype(X0), typeof(T), typeof(rho))

    fvec = zeros(T_out, 6)
    orders = X0[1:5]
    mu_B = X0[6]
    mus = [1/3*mu_B, 1/3*mu_B, 1/3*mu_B]
    fvec[1:5] = dOmega_dorder(orders, mus, T, ints)
    rho_now = - sum(dOmega_dmus(orders, mus, T, ints)) / rho0
    fvec[6] = rho_now - rho
    return fvec
end


function Tmu(X0, mu_B, T, ints)
    fWrapper(Xs) = Quark_mu(Xs, mu_B, T, ints)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    NewX = res.zero
    return NewX
end




function Trho(X0, T, rho, ints)
    fWrapper(Xs) = Quark_rho(Xs, T, rho, ints)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    NewX = res.zero
    return NewX
end


function get_EOS(NewX, T, rho, theta, Gv, ints)
    orders = NewX[1:8]
    mu_B,mu_Q = NewX[9:10]
    mu_u = 1/3*mu_B + 2/3 * mu_Q
    mu_d = 1/3*mu_B - 1/3 * mu_Q
    mu_s = 1/3*mu_B - 1/3 * mu_Q
    mu_e = - mu_Q
    mu_mu = - mu_Q
    mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
    rhoB = 3 * rho * rho0
    mu_u_eff = mu_u - 2 * Gv * rhoB
    mu_d_eff = mu_d - 2 * Gv * rhoB
    mu_s_eff = mu_s - 2 * Gv * rhoB
    P = - Omega(orders, mus, T, theta, rho, Gv, ints)
    rho_u, rho_d, rho_s, rho_e, rho_mu = -dOmgea_dmus(orders, mus, T, theta, rho, Gv, ints)
    S = - dOmgea_dT(orders, mus, T, theta, rho, Gv, ints)
    rho_mu_eff = rho_u * mu_u_eff + rho_d * mu_d_eff + rho_s * mu_s_eff + rho_e * mu_e + rho_mu * mu_mu
    E = - P + T * S + rho_mu_eff
    return [rho, P, E]
end