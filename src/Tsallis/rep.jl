using Revise

includet("pnjl_TT.jl")
using CSV 
using DataFrames
using FiniteDifferences
using ForwardDiff

function SolveOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    X0_typed = convert.(promote_type(eltype(X0), typeof(T), typeof(mu_B)), X0)
    
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B, eB, q, gauss_nodes, zeta_nodes)

    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    NewX = res.zero
    phi = NewX[1:3]
    Phi1 = NewX[4]
    Phi2 = NewX[5]
    return -Omega(phi, Phi1, Phi2, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
end


function DmuOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    return ForwardDiff.derivative(x -> SolveOmega(X0, T, x, eB, q, gauss_nodes, zeta_nodes), mu_B)
end

function Dmu2Omega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    # d^2 P / dmu^2
    return ForwardDiff.derivative(x -> DmuOmega(X0, T, x, eB, q, gauss_nodes, zeta_nodes), mu_B)
end

function Dmu3Omega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    # d^3 P / dmu^3
    return ForwardDiff.derivative(x -> Dmu2Omega(X0, T, x, eB, q, gauss_nodes, zeta_nodes), mu_B)
end

function Dmu4Omega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    # d^4 P / dmu^4
    return ForwardDiff.derivative(x -> Dmu3Omega(X0, T, x, eB, q, gauss_nodes, zeta_nodes), mu_B)
end


function DTOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    # dP / dT
    return ForwardDiff.derivative(x -> SolveOmega(X0, x, mu_B, eB, q, gauss_nodes, zeta_nodes), T)
end

function DTTOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    # d^2 P / dT^2
    return ForwardDiff.derivative(x -> DTOmega(X0, x, mu_B, eB, q, gauss_nodes, zeta_nodes), T)
end

function DmuTOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    # d^2 P / dmu dT
    return ForwardDiff.derivative(x -> DmuOmega(X0, x, mu_B, eB, q, gauss_nodes, zeta_nodes), T)
end







function Ther_Rep(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes, P0)

    P = SolveOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes) - P0
    chi_T = DTOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    chi_mu = DmuOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    chi_mumu = Dmu2Omega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    chi_TT = DTTOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    chi_muT = DmuTOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)

    E = -P + T*chi_T + mu_B*chi_mu
    

    v_n_2 = (chi_T * chi_mumu - chi_mu * chi_muT) / (T * (chi_mumu * chi_TT - chi_muT^2))
    v_s_2 = (chi_T * chi_muT - chi_mu * chi_TT) / (mu_B * (chi_muT^2 - chi_TT * chi_mumu))
    E = T * chi_T + mu_B * chi_mu - P
    v_2 = (v_n_2 * T * chi_T + v_s_2 * mu_B * chi_mu) / (P + E)


    return [T*197.33, mu_B*197.33, P, E, v_2]
end


# 涨落
function Fluctuations(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes)

    chi_mu = DmuOmega(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes) / T^3     # n = dP/dmu
    chi2_mu2 = Dmu2Omega(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes) / T^2 # d2P/dmu2
    chi3_mu3 = Dmu3Omega(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes) /T # d3P/dmu3
    chi4_mu4 = Dmu4Omega(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes) # d4P/dmu4
    P = SolveOmega(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    chi21 = chi2_mu2 / chi_mu
    chi31 = chi3_mu3 / chi_mu
    chi42 = chi4_mu4 / chi2_mu2

    return [T*197.33, mu_B*197.33, P, chi_mu, chi2_mu2, chi3_mu3, chi4_mu4]
end

function NumDmuOmega(X0, T, mu_B, eB, q, gauss_nodes, zeta_nodes)

    method1 = central_fdm(5, 1)  # 4阶精度
    method2 = central_fdm(7, 2)  # 4阶精度
    method3 = central_fdm(7, 3)  # 4阶精度
    method4 = central_fdm(9, 4)  # 4阶精度
    f(x) = SolveOmega(X0, T, x, eB, q, gauss_nodes, zeta_nodes)
    res1 = method1(f, mu_B)
    res2 = method2(f, mu_B)
    res3 = method3(f, mu_B)
    res4 = method4(f, mu_B)
    return [res1, res2, res3, res4]
end




function num_Flu(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    P = SolveOmega(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    res = NumDmuOmega(NewX, T, mu_B, eB, q, gauss_nodes, zeta_nodes)
    chi_mu = res[1] / T^3     # n = dP/dmu
    chi2_mu2 = res[2] / T^2 # d2P/dmu2
    chi3_mu3 = res[3] /T # d3P/dmu3
    chi4_mu4 = res[4] # d4P/dmu4
    return [T*197.33, mu_B*197.33, P, chi_mu, chi2_mu2, chi3_mu3, chi4_mu4]
end