includet("Pnjl_TT.jl")
using CSV 
using DataFrames

function SolveOmega(X0, T, mu_B, q, ints)
    X0_typed = convert.(promote_type(eltype(X0), typeof(T), typeof(mu_B)), X0)
    
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B, q, ints)

    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    NewX = res.zero
    phi = NewX[1:3]
    Phi1 = NewX[4]
    Phi2 = NewX[5]
    return -Omega(phi, Phi1, Phi2, T, mu_B, q, ints)
end


function DmuOmega(X0, T, mu_B, q, ints)
    return ForwardDiff.derivative(x -> SolveOmega(X0, T, x, q, ints), mu_B)
end

function Dmu2Omega(X0, T, mu_B, q, ints)
    # d^2 P / dmu^2
    return ForwardDiff.derivative(x -> DmuOmega(X0, T, x, q, ints), mu_B)
end

function Dmu3Omega(X0, T, mu_B, q, ints)
    # d^3 P / dmu^3
    return ForwardDiff.derivative(x -> Dmu2Omega(X0, T, x, q, ints), mu_B)
end

function Dmu4Omega(X0, T, mu_B, q, ints)
    # d^4 P / dmu^4
    return ForwardDiff.derivative(x -> Dmu3Omega(X0, T, x, q, ints), mu_B)
end

function  DTOmega(X0, T, mu_B, q, ints)
    # dP / dT
    return ForwardDiff.derivative(x -> SolveOmega(X0, x, mu_B, q, ints), T)
end

function DTTOmega(X0, T, mu_B, q, ints)
    # d^2 P / dT^2
    return ForwardDiff.derivative(x -> DTOmega(X0, x, mu_B, q, ints), T)
end

function DmuTOmega(X0, T, mu_B, q, ints)
    # d^2 P / dmu dT
    return ForwardDiff.derivative(x -> DmuOmega(X0, x, mu_B, q, ints), T)
end






# 涨落
function Fluctuations(NewX, T, mu_B, P0, q, ints)

    chi_mu = DmuOmega(NewX, T, mu_B, q, ints) / T^3     # n = dP/dmu
    chi2_mu2 = Dmu2Omega(NewX, T, mu_B, q, ints) / T^2 # d2P/dmu2
    chi3_mu3 = Dmu3Omega(NewX, T, mu_B, q, ints) /T # d3P/dmu3
    chi4_mu4 = Dmu4Omega(NewX, T, mu_B, q, ints) # d4P/dmu4
    P = SolveOmega(NewX, T, mu_B, q, ints) - P0
    chi21 = chi2_mu2 / chi_mu
    chi31 = chi3_mu3 / chi_mu
    chi42 = chi4_mu4 / chi2_mu2

    return [T*197.33, mu_B*197.33, P, chi_mu, chi2_mu2, chi3_mu3, chi4_mu4, mu_B/T, P/T^4]
end



function Ther_Rep(X0, T, mu_B,  P0, q, ints)

    P = SolveOmega(X0, T, mu_B, q, ints) - P0
    chi_T = DTOmega(X0, T, mu_B, q, ints)
    S = chi_T
    chi_mu = DmuOmega(X0, T, mu_B, q, ints)
    rho_B = chi_mu
    chi_mumu = Dmu2Omega(X0, T, mu_B, q, ints)
    chi_TT = DTTOmega(X0, T, mu_B, q, ints)
    chi_muT = DmuTOmega(X0, T, mu_B, q, ints)

    E = -P + T*chi_T + mu_B*chi_mu
    # 热容


    CV = T * (chi_TT - chi_muT^2 / chi_mumu)
    

    #CP = T * (chi_TT - 2*chi_T*chi_muT/chi_mu + (chi_T/chi_mu)^2 * chi_mumu)
    CP = T * (chi_TT - (chi_T * chi_muT) / chi_mu)
    TA = E-3*P


    v_n_2 = (chi_T * chi_mumu - chi_mu * chi_muT) / (T * (chi_mumu * chi_TT - chi_muT^2)) # 等密声速平方
    v_s_2 = (chi_T * chi_muT - chi_mu * chi_TT) / (mu_B * (chi_muT^2 - chi_TT * chi_mumu)) # 等熵声速平方
    # 等 s/rho
    v_2 = (v_n_2 * T * chi_T + v_s_2 * mu_B * chi_mu) / (P + E) # P+E = T*chi_T + mu_B*chi_mu
    
    if !isfinite(v_2)
        v_2 = 0.0
    end

    return [T*197.33, mu_B*197.33, P/T^4, E/T^4, TA/T^4, S/T^3, CP/T^3, CV/T^3, v_2]
    #return [T*197.33, mu_B*197.33, P*hc, E*hc, TA*hc, S, rho_B, CV, v_2]
end

