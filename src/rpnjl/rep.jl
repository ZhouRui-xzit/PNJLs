include("../../src/rpnjl/rpnjl.jl")
using CSV 
using DataFrames


function SolveOmega(X0, mu_B, T, ints)
    X0_typed = convert.(promote_type(eltype(X0), typeof(T), typeof(mu_B)), X0)

    fWrapper(Xs) = Quark_mu(Xs, mu_B, T, ints)

    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    NewX = res.zero
    orders = NewX[1:5]
    mus = [1/3*mu_B, 1/3*mu_B, 1/3*mu_B]

    return -Omega(orders, mus, T, ints)
end



function DmuOmega(X0, mu_B, T, ints)
    return ForwardDiff.derivative(x -> SolveOmega(X0, x, T, ints), mu_B)
end

function Dmu2Omega(X0, mu_B, T, ints)
    # d^2 P / dmu^2
    return ForwardDiff.derivative(x -> DmuOmega(X0, x, T, ints), mu_B)
end

function Dmu3Omega(X0, mu_B, T, ints)
    # d^3 P / dmu^3
    return ForwardDiff.derivative(x -> Dmu2Omega(X0, x, T, ints), mu_B)
end

function Dmu4Omega(X0, mu_B, T, ints)
    # d^4 P / dmu^4
    return ForwardDiff.derivative(x -> Dmu3Omega(X0, x, T, ints), mu_B)
end


# 涨落
function Fluctuations(NewX, mu_B, T, ints)

    chi_mu = DmuOmega(NewX, mu_B, T, ints) / T^3
    chi2_mu2 = Dmu2Omega(NewX, mu_B, T, ints) / T^2
    chi3_mu3 = Dmu3Omega(NewX, mu_B, T, ints) / T 
    chi4_mu4 = Dmu4Omega(NewX, mu_B, T, ints) 
    P = (SolveOmega(NewX, mu_B, T, ints) - 26.78879) / T^4

    chi21 = chi2_mu2 / chi_mu
    chi32 = chi3_mu3 / chi2_mu2
    chi42 = chi4_mu4 / chi2_mu2



    return [T*197.33, mu_B*197.33, P, chi_mu, chi2_mu2, chi3_mu3, chi4_mu4]
end