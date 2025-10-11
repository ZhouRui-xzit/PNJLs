include("pnjl_magnetic.jl")
using CSV 
using DataFrames
using FiniteDifferences


function SolveOmega(X0, T, mu_B, eB, ints)
    X0_typed = convert.(promote_type(eltype(X0), typeof(T), typeof(mu_B)), X0)
    
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B, eB, ints)

    res = nlsolve(fWrapper, X0_typed, autodiff=:finite)
    NewX = res.zero
    phi = NewX[1:3]
    Phi1 = NewX[4]
    Phi2 = NewX[5]
    return -Omega(phi, Phi1, Phi2, T, mu_B, eB, ints)
end


function DmuOmega(X0, T, mu_B, eB, ints)
    return ForwardDiff.derivative(x -> SolveOmega(X0, T, x, eB, ints), mu_B)
end

function Dmu2Omega(X0, T, mu_B, eB, ints)
    # d^2 P / dmu^2
    return ForwardDiff.derivative(x -> DmuOmega(X0, T, x, eB, ints), mu_B)
end

function Dmu3Omega(X0, T, mu_B, eB, ints)
    # d^3 P / dmu^3
    return ForwardDiff.derivative(x -> Dmu2Omega(X0, T, x, eB, ints), mu_B)
end

function Dmu4Omega(X0, T, mu_B, eB, ints)
    # d^4 P / dmu^4
    return ForwardDiff.derivative(x -> Dmu3Omega(X0, T, x, eB, ints), mu_B)
end


function DTOmega(X0, T, mu_B, eB, ints)
    # dP / dT
    return ForwardDiff.derivative(x -> SolveOmega(X0, x, mu_B, eB, ints), T)
end

function DTTOmega(X0, T, mu_B, eB, ints)
    # d^2 P / dT^2
    return ForwardDiff.derivative(x -> DTOmega(X0, x, mu_B, eB, ints), T)
end

function DmuTOmega(X0, T, mu_B, eB, ints)
    # d^2 P / dmu dT
    return ForwardDiff.derivative(x -> DmuOmega(X0, x, mu_B, eB, ints), T)
end

function numDerv(X0, T, mu_B, eB, q)
    method1 = FiniteDifferences.central_fdm(5, 1)
    method2 = FiniteDifferences.central_fdm(5, 2)
    method3 = FiniteDifferences.central_fdm(7, 3)
    method4 = FiniteDifferences.central_fdm(9, 4)
    f(x) = SolveOmega(X0, T, x, eB, q)
    chi_mu = method1(f, mu_B)
    chi2_mu2 = method2(f, mu_B)
    chi3_mu3 = method3(f, mu_B)
    chi4_mu4 = method4(f, mu_B)
    return chi_mu, chi2_mu2, chi3_mu3, chi4_mu4
end




# 涨落
function Fluctuations(NewX, T, mu_B, eB, ints)

    chi_mu = DmuOmega(NewX, T, mu_B, eB, ints)     # n = dP/dmu
    chi2_mu2 = Dmu2Omega(NewX, T, mu_B, eB, ints) # d2P/dmu2
    chi3_mu3 = Dmu3Omega(NewX, T, mu_B, eB, ints) # d3P/dmu3
    chi4_mu4 = Dmu4Omega(NewX, T, mu_B, eB, ints) # d4P/dmu4

    chi21 = chi2_mu2 / chi_mu
    chi31 = chi3_mu3 / chi_mu
    chi42 = chi4_mu4 / chi2_mu2

    return [T*197.33, mu_B*197.33, chi21, chi31, chi42]
end


function Ther_Rep(X0, T, mu_B, eB, ints, P0)

    P = SolveOmega(X0, T, mu_B, eB, ints) - P0
    chi_T = DTOmega(X0, T, mu_B, eB, ints)
    chi_mu = DmuOmega(X0, T, mu_B, eB, ints)
    chi_mumu = Dmu2Omega(X0, T, mu_B, eB, ints)
    chi_TT = DTTOmega(X0, T, mu_B, eB, ints)
    chi_muT = DmuTOmega(X0, T, mu_B, eB, ints)

    E = -P + T*chi_T + mu_B*chi_mu
    #Tr_A = (E - 3*P) /T^4 

    v_n_2 = (chi_T * chi_mumu - chi_mu * chi_muT) / (T * (chi_mumu * chi_TT - chi_muT^2))
    v_s_2 = (chi_T * chi_muT - chi_mu * chi_TT) / (mu_B * (chi_muT^2 - chi_TT * chi_mumu))
    E = T * chi_T + mu_B * chi_mu - P
    v_2 = (v_n_2 * T * chi_T + v_s_2 * mu_B * chi_mu) / (P + E)


    return [T*197.33, mu_B*197.33, P, E, v_n_2, v_s_2, v_2]
end


function numF(NewX, T, mu_B, eB, ints)
    chi_mu, chi2_mu2, chi3_mu3, chi4_mu4 = numDerv(NewX, T, mu_B, eB, ints)
    chi_mu = chi_mu / T^3 
    chi2_mu2 = chi2_mu2 / T^2
    chi3_mu3 = chi3_mu3 / T
    chi4_mu4 = chi4_mu4
    chi21 = chi2_mu2 / chi_mu
    chi31 = chi3_mu3 / chi_mu
    chi42 = chi4_mu4 / chi2_mu2
    return [T*197.33, mu_B*197.33, chi21, chi31, chi42]
end
