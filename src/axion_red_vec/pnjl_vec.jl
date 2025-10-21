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
    p1, w1 = gauleg(0.0,Lambda_f, p_num)
    p2, w2 = gauleg(0.0,20.0, p_num)
    w1 = w1 .* p1.^2 ./ (2*pi^2)
    w2 = w2 .* p2.^2 ./ (2*pi^2)
    int1 = (p1, w1)
    int2 = (p2, w2)
    return int1, int2
end





function Omega(orders, mus, T, theta, rho, Gv, ints)

    p1, w1 = ints[1]
    p2, w2 = ints[2]

    phi = orders[1:6]
    Phi1 = orders[7]
    Phi2 = orders[8]

    rhoB = 3 * rho * rho0

    chi = chiral(phi, theta)
    U = calc_U(T, Phi1, Phi2)
    Omega_q = 0.0
    for flavor = 1:3
        mu = mus[flavor] - 2 * Gv * rhoB
        mass = Mass(phi, theta)[flavor]
        vacuum_term =  calculate_vacuum_term(p1, w1, mass)
        thermal_term = calculate_thermal_term(p2, w2, mass, T, mu, Phi1, Phi2)
        Omega_q += vacuum_term + thermal_term
    end

    mu_e = mus[4]
    mu_mu = mus[5]
    E_e = sqrt.(p2.^2 .+ m_e^2)
    E_mu = sqrt.(p2.^2 .+ m_mu^2)

    Omega_e = -2*T * sum(w2 .* fermion.(E_e .- mu_e, T) .+ w2 .* fermion.(E_e .+ mu_e, T))
    Omega_mu = -2*T * sum(w2 .* fermion.(E_mu .- mu_mu, T) .+ w2 .* fermion.(E_mu .+ mu_mu, T))


    Omega_total = chi + U + Omega_q - Gv * rhoB^2 + Omega_e + Omega_mu
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
function dOmgea_dorder(orders, mus, T, theta, rho, Gv, ints)
    return ForwardDiff.gradient(x -> Omega(x, mus, T, theta, rho, Gv, ints), orders)
end
 


# 对mu求导
function dOmgea_dmus(orders, mus, T, theta, rho, Gv, ints)
    return ForwardDiff.gradient(x -> Omega(orders, x, T, theta, rho, Gv, ints), mus)
end

# 对T求导 
function dOmgea_dT(orders, mus, T, theta, rho, Gv, ints)
    return ForwardDiff.derivative(x -> Omega(orders, mus, x, theta, rho, Gv, ints), T)
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


function fphi(x, T, Phi1, Phi2)
    term1 = Phi1 * exp.(-x ./ T) + 2 * Phi2 * exp.(-2.0 * x ./ T) +exp.(-3.0 * x ./ T)
    term2 = 1.0 + 3.0 * Phi1 * exp.(-x ./ T) + 3.0 * Phi2 * exp.(-2.0 * x ./ T) + exp.(-3.0 * x ./ T)
    over = -3.0 * x / T 

    if over > 307.0
        return 1.0
    elseif over < -307.0
        return 0.0
    else 
        return term1 ./ term2
    end

end

function fphibar(x, T, Phi1, Phi2)
    term1 = Phi2 * exp.(-x ./ T) + 2 * Phi1 * exp.(-2.0 * x ./ T) + exp.(-3.0 * x ./ T)
    term2 = 1.0 + 3.0 * Phi2 * exp.(-x ./ T) + 3.0 * Phi1 * exp.(-2.0 * x ./ T) + exp.(-3.0 * x ./ T)

    over = -3.0 * x / T

    if over > 307.0
        return 1.0
    elseif over < -307.0
        return 0.0
    else 
        return term1 ./ term2
    end

end






function fermion(x, T)
    term1 = log.(1 .+ exp.(-x ./ T))
    return term1
end

function n_fermion(x, T)
    over = x / T
    if over > 307.0
        return 0.0
    else 
        return 1.0 / (1.0 + exp.(over))
    end
end



function chiral(phi, theta)
    @inbounds begin
        sigma = phi[1:3]
        eta   = phi[4:6]
        sigma_u, sigma_d, sigma_s = sigma
        eta_u,   eta_d,   eta_s   = eta
    end

    # 二次项
    term1 = 2 * G_f * (sum(sigma .^ 2) + sum(eta .^ 2))

    # 't Hooft 六次项（含 theta 耦合）
    term2 = 4 * K_f * (
          cos(theta) * (sigma_u * sigma_d * sigma_s)
        + sin(theta) * (eta_u   * eta_d   * eta_s)
        - cos(theta) * (sigma_u * eta_d   * eta_s + eta_u * eta_d * sigma_s + eta_u * eta_s * sigma_d)
        - sin(theta) * (sigma_u * sigma_s * eta_d + sigma_u * sigma_d * eta_s + sigma_d * sigma_s * eta_u)
    )

    return term1 + term2
end




function Mass(phi, theta)
    """计算三种夸克的有效质量"""
    sigma = phi[1:3]
    eta   = phi[4:6]

    mu_s = m0[1] + 4 * G_f * sigma[1] + 2 * K_f * (
            cos(theta) * (sigma[2] * sigma[3] - eta[2] * eta[3]) -
            sin(theta) * (sigma[2] * eta[3] + eta[2] * sigma[3])
    )

    mu_ps = 4 * G_f * eta[1] - 2 * K_f * (
            cos(theta) * (sigma[2] * eta[3] + eta[2] * sigma[3]) +
            sin(theta) * (sigma[2] * sigma[3] - eta[2] * eta[3])
    )
    
    mass_u = sqrt(mu_s^2 + mu_ps^2)

    md_s = m0[2] + 4 * G_f * sigma[2] + 2 * K_f * (
            cos(theta) * (sigma[3] * sigma[1] - eta[3] * eta[1]) -
            sin(theta) * (sigma[3] * eta[1] + eta[3] * sigma[1])
    )
    
    md_ps = 4 * G_f * eta[2] - 2 * K_f * (
            cos(theta) * (sigma[3] * eta[1] + eta[3] * sigma[1]) +
            sin(theta) * (sigma[3] * sigma[1] - eta[3] * eta[1])
    )
    
    mass_d = sqrt(md_s^2 + md_ps^2)
    
    ms_s = m0[3] + 4 * G_f * sigma[3] + 2 * K_f * (
            cos(theta) * (sigma[1] * sigma[2] - eta[1] * eta[2]) -
            sin(theta) * (sigma[1] * eta[2] + eta[1] * sigma[2])
    )
    
    ms_ps = 4 * G_f * eta[3] - 2 * K_f * (
            cos(theta) * (sigma[1] * eta[2] + eta[1] * sigma[2]) +
            sin(theta) * (sigma[1] * sigma[2] - eta[1] * eta[2])
    )
    
    mass_s = sqrt(ms_s^2 + ms_ps^2)
    

    return [mass_u, mass_d, mass_s]
end

E_of_p(p, mass) = sqrt(p*p + mass*mass)

function calc_U(T, Phi1, Phi2)
    # 对数型 Polyakov 势
    Ta = a0 + a1 * (T0/T) + a2 * (T0/T)^2
    Tb = b3 * (T0 / T)^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    U = T^4 * (-0.5 * Ta * Phi2 * Phi1 + Tb * NaNMath.log(value))
    return U
end


function rho_mu_eff(mu_B, mu_Q, T,theta, NewX, rho, Gv, ints2)
    mu_u = 1/3*mu_B + 2/3 * mu_Q
    mu_d = 1/3*mu_B - 1/3 * mu_Q
    mu_s = 1/3*mu_B - 1/3 * mu_Q
    rho = 3 * rho * rho0  

    mu_u_eff = mu_u - 2 * Gv * rho
    mu_d_eff = mu_d - 2 * Gv * rho
    mu_s_eff = mu_s - 2 * Gv * rho
    mues = [mu_u_eff, mu_d_eff, mu_s_eff]
    mu_e = - mu_Q
    mu_mu = - mu_Q
    

    phi = NewX[1:6]
    Phi1, Phi2 = NewX[7:8]
    p, w = ints2
    Masses = Mass(phi, theta)

    rho_eff = 0.0
    rho_flavor = [0.0, 0.0, 0.0]
  
    for flavor = 1:3
        E = sqrt.(p.^2 .+ Masses[flavor]^2)
        E_minus = E .- mues[flavor]
        E_plus = E .+ mues[flavor]
        f_dist = fphi.(E_minus, T, Phi1, Phi2)
  
        fbar_dist = fphibar.(E_plus, T, Phi1, Phi2)
        rho_flavor[flavor] = sum(6 * w .* (f_dist .- fbar_dist))
        rho_eff += rho_flavor[flavor] * mues[flavor]
    end

    E_e = sqrt.(p.^2 .+ m_e^2)
    E_mu = sqrt.(p.^2 .+ m_mu^2)
    rho_e = sum(2 * w .* n_fermion.(E_e .- mu_e, T) .- 2 * w .* n_fermion.(E_e .+ mu_e, T))
    rho_mu = sum(2 * w .* n_fermion.(E_mu .- mu_mu, T) .- 2 * w .* n_fermion.(E_mu .+ mu_mu, T))
    rho_eff += rho_e * mu_e + rho_mu * mu_mu
    return rho_eff, rho_flavor[1], rho_flavor[2], rho_flavor[3], rho_e, rho_mu
end









"""
    计算给定T和mu_B下的序参量
    X0 = [phi_u, phi_d, phi_s...., Phi1, Phi2, mu_Q, rho] 初始猜测值
"""

function Quark_mu(X0, mu_B, T, theta, Gv, ints)

    T_out = promote_type(eltype(X0), typeof(T), typeof(mu_B))
    orders = X0[1:8]
    rho = X0[9]
    mu_Q = X0[10]

    mu_u = 1/3*mu_B + 2/3 * mu_Q
    mu_d = 1/3*mu_B - 1/3 * mu_Q
    mu_s = 1/3*mu_B - 1/3 * mu_Q
    mu_e = - mu_Q
    mu_mu = - mu_Q
    mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
    fvec = zeros(T_out, 10)
    fvec[1:8] = dOmgea_dorder(orders, mus, T, theta, rho, Gv, ints)
    rho_u, rho_d, rho_s, rho_e, rho_mu = -dOmgea_dmus(orders, mus, T, theta, rho, Gv, ints)
    fvec[9] = rho_u + rho_d + rho_s - 3 * rho * rho0  # 体积密度转换为数密度
    fvec[10] = (2/3) * rho_u - (1/3) * rho_d - (1/3) * rho_s - rho_e - rho_mu  # 电荷守恒
    return fvec
end




"""
    计算给定T和rho_B下的序参量
    X0 = [phi_u, phi_d, phi_s, Phi1, Phi2, muB, muQ] 初始猜测值
"""

function Quark_rho(X0, T, rho, theta, Gv, ints)
    T_out = promote_type(eltype(X0), typeof(T), typeof(rho))

    fvec = zeros(T_out, 10)
    orders = X0[1:8]
    muB = X0[9]
    muQ = X0[10]
    mu_u = 1/3*muB + 2/3 * muQ
    mu_d = 1/3*muB - 1/3 * muQ
    mu_s = 1/3*muB - 1/3 * muQ
    mu_e = - muQ
    mu_mu = - muQ
    mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]
    fvec[1:8] = dOmgea_dorder(orders, mus, T, theta, rho, Gv, ints)
    rho_u, rho_d, rho_s, rho_e, rho_mu = -dOmgea_dmus(orders, mus, T, theta, rho, Gv, ints)
    fvec[9] = rho_u + rho_d + rho_s - 3 * rho * rho0  # 体积密度转换为数密度
    fvec[10] = (2/3) * rho_u - (1/3) * rho_d - (1/3) * rho_s - rho_e - rho_mu  # 电荷守恒

    return fvec
end


function Tmu(X0, mus, T, theta, Gv, ints)
    fWrapper(Xs) = Quark_mu(Xs, mus, T, theta, Gv, ints)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    NewX = res.zero
    return NewX
end




function Trho(X0, T, rho, theta, Gv, ints)
    fWrapper(Xs) = Quark_rho(Xs, T, rho, theta, Gv, ints)
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