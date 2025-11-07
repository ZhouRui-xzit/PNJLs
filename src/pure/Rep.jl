includet("Pnjl_pure.jl")
using CSV 
using DataFrames

function SolveOmega(X0, T, mu_B,)
    X0_typed = convert.(promote_type(eltype(X0), typeof(T), typeof(mu_B)), X0)
    
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B)

    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    NewX = res.zero
    phi = NewX[1:3]
    Phi1 = NewX[4]
    Phi2 = NewX[5]
    return -Omega(phi, Phi1, Phi2, T, mu_B)
end


function DmuOmega(X0, T, mu_B)
    return ForwardDiff.derivative(x -> SolveOmega(X0, T, x), mu_B)
end

function Dmu2Omega(X0, T, mu_B)
    # d^2 P / dmu^2
    return ForwardDiff.derivative(x -> DmuOmega(X0, T, x), mu_B)
end

function Dmu3Omega(X0, T, mu_B)
    # d^3 P / dmu^3
    return ForwardDiff.derivative(x -> Dmu2Omega(X0, T, x), mu_B)
end

function Dmu4Omega(X0, T, mu_B)
    # d^4 P / dmu^4
    return ForwardDiff.derivative(x -> Dmu3Omega(X0, T, x), mu_B)
end

# 涨落
function Fluctuations(NewX, T, mu_B)

    chi_mu = DmuOmega(NewX, T, mu_B) / T^3     # n = dP/dmu
    chi2_mu2 = Dmu2Omega(NewX, T, mu_B) / T^2 # d2P/dmu2
    chi3_mu3 = Dmu3Omega(NewX, T, mu_B) /T # d3P/dmu3
    chi4_mu4 = Dmu4Omega(NewX, T, mu_B) # d4P/dmu4
    P = SolveOmega(NewX, T, mu_B) / T^4 
    chi21 = chi2_mu2 / chi_mu
    chi31 = chi3_mu3 / chi_mu
    chi42 = chi4_mu4 / chi2_mu2

    return [T*197.33, mu_B*197.33, P, chi_mu, chi2_mu2, chi3_mu3, chi4_mu4]
end


