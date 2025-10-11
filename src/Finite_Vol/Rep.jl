include("pnjl_FV.jl")
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

function Solve_R(X0, beta, gamma, ints)
    X0_typed = convert.(promote_type(eltype(X0), typeof(beta), typeof(gamma)), X0)
    
    fWrapper(Xs) = Quark_mu(Xs, 1/beta, -gamma/beta, ints)

    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    NewX = res.zero
    phi = NewX[1:3]
    Phi1 = NewX[4]
    Phi2 = NewX[5]
    return -Omega(phi, Phi1, Phi2, 1/beta, -gamma/beta, ints) * beta 
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

# 涨落
function Fluctuations(NewX, T, mu_B, ints)

    chi_mu = DmuOmega(NewX, T, mu_B, ints) / T^3     # n = dP/dmu
    chi2_mu2 = Dmu2Omega(NewX, T, mu_B, ints) / T^2 # d2P/dmu2
    chi3_mu3 = Dmu3Omega(NewX, T, mu_B, ints) /T # d3P/dmu3
    chi4_mu4 = Dmu4Omega(NewX, T, mu_B, ints) # d4P/dmu4

    chi21 = chi2_mu2 / chi_mu
    chi31 = chi3_mu3 / chi_mu
    chi42 = chi4_mu4 / chi2_mu2

    return [T*197.33, mu_B*197.33, chi21, chi31, chi42]
end


function Riemann(X0, beta ,gamma, ints) 
    # beta = 1/T, gamma = -mu_B/T

end