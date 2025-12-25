includet("Pnjl_pure.jl")
using CSV 
using DataFrames

function SolveOmega(X0, T, mu_B, ints)
    X0_typed = convert.(promote_type(eltype(X0), typeof(T), typeof(mu_B)), X0)
    
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B, ints)

    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    NewX = res.zero
    phi = NewX[1:3]
    Phi1 = NewX[4]
    Phi2 = NewX[5]
    return -Omega(phi, Phi1, Phi2, T, mu_B, ints)
end


function DmuOmega(X0, T, mu_B, ints)
    return ForwardDiff.derivative(x -> SolveOmega(X0, T, x, ints), mu_B)
end

function Dmu2Omega(X0, T, mu_B, ints)
    # d^2 P / dmu^2
    return ForwardDiff.derivative(x -> DmuOmega(X0, T, x, ints), mu_B)
end

function Dmu3Omega(X0, T, mu_B, ints)
    # d^3 P / dmu^3
    return ForwardDiff.derivative(x -> Dmu2Omega(X0, T, x, ints), mu_B)
end

function Dmu4Omega(X0, T, mu_B, ints)
    # d^4 P / dmu^4
    return ForwardDiff.derivative(x -> Dmu3Omega(X0, T, x, ints), mu_B)
end

function  DTOmega(X0, T, mu_B, ints)
    # dP / dT
    return ForwardDiff.derivative(x -> SolveOmega(X0, x, mu_B, ints), T)
end

function DTTOmega(X0, T, mu_B, ints)
    # d^2 P / dT^2
    return ForwardDiff.derivative(x -> DTOmega(X0, x, mu_B, ints), T)
end

function DmuTOmega(X0, T, mu_B, ints)
    # d^2 P / dmu dT
    return ForwardDiff.derivative(x -> DmuOmega(X0, x, mu_B, ints), T)
end






# 涨落
function Fluctuations(NewX, T, mu_B, ints)

    chi_mu = DmuOmega(NewX, T, mu_B, ints) / T^3     # n = dP/dmu
    chi2_mu2 = Dmu2Omega(NewX, T, mu_B, ints) / T^2 # d2P/dmu2
    chi3_mu3 = Dmu3Omega(NewX, T, mu_B, ints) /T # d3P/dmu3
    chi4_mu4 = Dmu4Omega(NewX, T, mu_B, ints) # d4P/dmu4
    P = SolveOmega(NewX, T, mu_B, ints)
    chi21 = chi2_mu2 / chi_mu
    chi31 = chi3_mu3 / chi_mu
    chi42 = chi4_mu4 / chi2_mu2

    return [T*197.33, mu_B*197.33, P, chi_mu, chi2_mu2, chi3_mu3, chi4_mu4]
end



function Ther_Rep(X0, T, mu_B,  P0, ints)

    P = SolveOmega(X0, T, mu_B, ints) - P0
    chi_T = DTOmega(X0, T, mu_B, ints)
    S = chi_T
    chi_mu = DmuOmega(X0, T, mu_B, ints)
    rho_B = chi_mu
    chi_mumu = Dmu2Omega(X0, T, mu_B, ints)
    chi_TT = DTTOmega(X0, T, mu_B, ints)
    chi_muT = DmuTOmega(X0, T, mu_B, ints)

    E = -P + T*chi_T + mu_B*chi_mu
    # 热容


    CV = T * (chi_TT - chi_muT^2 / chi_mumu)
    

    #CP = T * (chi_TT - 2*chi_T*chi_muT/chi_mu + (chi_T/chi_mu)^2 * chi_mumu)
    TA = E-3*P


    v_n_2 = (chi_T * chi_mumu - chi_mu * chi_muT) / (T * (chi_mumu * chi_TT - chi_muT^2)) # 等密声速平方
    v_s_2 = (chi_T * chi_muT - chi_mu * chi_TT) / (mu_B * (chi_muT^2 - chi_TT * chi_mumu)) # 等熵声速平方
    # 等 s/rho
    v_2 = (v_n_2 * T * chi_T + v_s_2 * mu_B * chi_mu) / (P + E) # P+E = T*chi_T + mu_B*chi_mu
    
    if !isfinite(v_2)
        v_2 = 0.0
    end

    return [T*197.33, mu_B*197.33, P/T^4, E/T^4, TA/T^4, S/T^3, rho_B, CV/T^3, v_2]
    #return [T*197.33, mu_B*197.33, P*hc, E*hc, TA*hc, S, rho_B, CV, v_2]
end




function alpha_S(T, mu)
    # 定义 Log 的参数 L，避免重复书写和计算，减少出错概率
    # 公式中是对数内的项： (T/Lambda_T) * sqrt(1 + (mu/(pi*T))^2)
    # 注意：T 和 mu 的单位必须与 Lambda_T 一致 (推荐都用 fm^-1)
    
    L_arg = (T / Lambda_T) * sqrt(1 + (mu / (pi * T))^2)
    L = log(L_arg)
    
    # 按照公式 (54)
    beta0 = 33 - 2 * Nf
    beta1 = 6 * (153 - 19 * Nf) # 注意系数: 3 * (...) / (...)^2，这里为了清晰拆开写
    
    term1 = (6 * pi) / (beta0 * L)
    
    # 公式第二部分: 1 - ...
    # 原公式中分子系数是 3*(153 - 19Nf)，分母是 (33 - 2Nf)^2
    factor = (3 * (153 - 19 * Nf)) / (beta0^2)
    
    term2 = 1 - factor * (log(2 * L) / L)
    
    return term1 * term2
end

function tau(T, mu)
    alphas = alpha_S(T, mu)
    
    # 按照公式 (53)
    # 分母 = 5.1 * T * alpha^2 * log(1/alpha) * (1 + 0.12*(2Nf + 1))
    denom = 5.1 * T * alphas^2 * log(1 / alphas) * (1 + 0.12 * (2 * Nf + 1))
    
    return 1 / denom
end


function quark_distribution(E, muq, T, Phi1, Phi2)
    n_f1 = 1 + 3 * Phi1 * exp(-(E - muq) / T) + 3 * Phi2 * exp(-2 * (E - muq) / T) + exp(-3 * (E - muq) / T)
    n_f2 = Phi1 * exp(-(E - muq) / T) + 2 * Phi2 * exp(-2 * (E - muq) / T) + exp(-3 * (E - muq) / T)
    return n_f2 / n_f1
end

function antiquark_distribution(E, muq, T, Phi1, Phi2)
    n_f1 = 1 + 3 * Phi2 * exp(-(E + muq) / T) + 3 * Phi1 * exp(-2 * (E + muq) / T) + exp(-3 * (E + muq) / T)
    n_f2 = Phi2 * exp(-(E + muq) / T) + 2 * Phi1 * exp(-2 * (E + muq) / T) + exp(-3 * (E + muq) / T)
    return n_f2 / n_f1
end




function trans_eff(X0, T, mu, ints)
    phi = X0[1:3]
    Phi1 = X0[4]
    Phi2 = X0[5]
    tau_val = tau(T, mu/3)
    masses = Mass(phi)
    eta = 0.0
    sigma_el =0.0
    muq = mu / 3.0
    p2, w2 = ints[2]
    e_quark = [2/3, -1/3, -1/3]  # 夸克电荷数组

    for flavor = 1:3

      
        E = sqrt.(p2.^2 .+ masses[flavor]^2)
        
        # 夸克分布函数 f⁰(1-f⁰)
        f_q = quark_distribution.(E, muq, T, Phi1, Phi2)
        f_qbar = antiquark_distribution.(E, muq, T, Phi1, Phi2)
        
        f_factor = f_q .* (1 .- f_q) .+ f_qbar .* (1 .- f_qbar)
        

        integrand_eta = p2.^4 ./ E.^2 .* f_factor .* w2  # 积分测度包含于w中

        integrand_sigma = e_quark[flavor]^2 * p2.^2 ./ E.^2 .* f_factor .* w2

        eta += sum(integrand_eta)
        sigma_el += sum(integrand_sigma)
    end
    eta *= 6*tau_val/(15*T)
    sigma_el *= 6*tau_val/(3*T)

    return [T*197.33, mu*197.33, eta, sigma_el/T]
end