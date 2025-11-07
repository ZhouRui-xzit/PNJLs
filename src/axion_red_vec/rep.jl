include("pnjl_vec.jl")
using CSV 
using DataFrames
using FiniteDifferences


function SolveOmega(NewX, mu_B, T, rho, theta, Gv, ints)
    # 确保 X0_typed 的类型与自动微分兼容
    X0_typed = convert.(promote_type(eltype(NewX), typeof(T), typeof(mu_B)), NewX)

    fWrapper(Xs) = Quark_rho(Xs, T, rho, theta, Gv, ints)
    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    orders = res.zero[1:8]
    mu_B = res.zero[9]
    mu_Q = res.zero[10]
    mu_u = 1/3*mu_B + 2/3 * mu_Q
    mu_d = 1/3*mu_B - 1/3 * mu_Q
    mu_s = 1/3*mu_B - 1/3 * mu_Q
    mu_e = - mu_Q
    mu_mu = - mu_Q
    mus = [mu_u, mu_d, mu_s, mu_e, mu_mu]

    return -Omega(orders, mus, T, theta, rho, Gv, ints)
end


function DmuOmega(NewX, mu_B, T, rho, theta, Gv, ints)
    return ForwardDiff.derivative(x -> SolveOmega(NewX, x, T, rho, theta, Gv, ints), mu_B)
end

function Dmu2Omega(NewX, mu_B, T, rho, theta, Gv, ints)
    # d^2 P / dmu^2
    return ForwardDiff.derivative(x -> DmuOmega(NewX, x, T, rho, theta, Gv, ints), mu_B)
end

function Dmu3Omega(NewX, mu_B, T, rho, theta, Gv, ints)
    # d^3 P / dmu^3
    return ForwardDiff.derivative(x -> Dmu2Omega(NewX, x, T, rho, theta, Gv, ints), mu_B)
end

function Dmu4Omega(NewX, mu_B, T, rho, theta, Gv, ints)
    # d^4 P / dmu^4
    return ForwardDiff.derivative(x -> Dmu3Omega(NewX, x, T, rho, theta, Gv, ints), mu_B)
end


function DTOmega(NewX, mu_B, T, rho, theta, Gv, ints)
    # dP / dT
    return ForwardDiff.derivative(x -> SolveOmega(NewX, mu_B, x, rho, theta, Gv, ints), T)
end

function DTTOmega(NewX, mu_B, T, rho, theta, Gv, ints)
    # d^2 P / dT^2
    return ForwardDiff.derivative(x -> DTOmega(NewX, mu_B, x, rho, theta, Gv, ints), T)
end

function DmuTOmega(NewX, mu_B, T, rho, theta, Gv, ints)
    # d^2 P / dmu dT
    return ForwardDiff.derivative(x -> DmuOmega(NewX, mu_B, x, rho, theta, Gv, ints), T)
end

function numDerv(X0, T, mu_B, q)
    method1 = FiniteDifferences.central_fdm(5, 1)
    method2 = FiniteDifferences.central_fdm(5, 2)
    method3 = FiniteDifferences.central_fdm(7, 3)
    method4 = FiniteDifferences.central_fdm(9, 4)
    f(x) = SolveOmega(X0, T, x, q)
    chi_mu = method1(f, mu_B)
    chi2_mu2 = method2(f, mu_B)
    chi3_mu3 = method3(f, mu_B)
    chi4_mu4 = method4(f, mu_B)
    return chi_mu, chi2_mu2, chi3_mu3, chi4_mu4
end



# 涨落
function Fluctuations(NewX, T, mu_B, q)

    chi_mu = DmuOmega(NewX, T, mu_B, q)     # n = dP/dmu
    chi2_mu2 = Dmu2Omega(NewX, T, mu_B, q) # d2P/dmu2
    chi3_mu3 = Dmu3Omega(NewX, T, mu_B, q) # d3P/dmu3
    chi4_mu4 = Dmu4Omega(NewX, T, mu_B, q) # d4P/dmu4
    P = SolveOmega(NewX, T, mu_B, q)
    chi21 = chi2_mu2 / chi_mu
    chi31 = chi3_mu3 / chi_mu
    chi42 = chi4_mu4 / chi2_mu2

    return [T*197.33, mu_B*197.33, P, chi_mu, chi2_mu2, chi3_mu3, chi4_mu4]
end



function Ther_Rep(NewX, mu_B, T, rho, theta, Gv, ints)
    




    chi_T =  DTOmega(NewX, mu_B, T, rho, theta, Gv, ints)
    chi_mu = DmuOmega(NewX, mu_B, T, rho, theta, Gv, ints)
    chi_TT = DTTOmega(NewX, mu_B, T, rho, theta, Gv, ints)
    chi_mumu = Dmu2Omega(NewX, mu_B, T, rho, theta, Gv, ints)
    chi_muT = DmuTOmega(NewX, mu_B, T, rho, theta, Gv, ints)


   #v2 = (chi_T * chi_mumu - chi_mu * chi_muT) / (T * (chi_mumu * chi_TT - chi_muT^2))
    #P = -SolveOmega(NewX, mu_B, T, rho, theta, Gv, ints)

    return [T*197.33, mu_B*197.33, chi_T, chi_mu, chi_TT, chi_mumu, chi_muT]
end
