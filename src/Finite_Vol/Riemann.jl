include("pnjl_FV.jl")
using CSV 
using DataFrames
using LinearAlgebra


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

function DbetaDR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> Solve_R(X0, x, gamma, ints), beta)
end

function Dbeta2DR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> DbetaDR(X0, x, gamma, ints), beta)
end

function Dbeta3DR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> Dbeta2DR(X0, x, gamma, ints), beta)
end


function DgammaDR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> Solve_R(X0, beta, x, ints), gamma)
end


function Dgamma2DR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> DgammaDR(X0, beta, x, ints), gamma)
end

function Dgamma3DR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> Dgamma2DR(X0, beta, x, ints), gamma)
end

function Dbeta_gammaDR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> DgammaDR(X0, x, gamma, ints), beta)
end

function Dbeta2_gammaDR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> Dbeta2DR(X0, beta, x, ints), gamma)
end

function Dgamma2_betaDR(X0, beta, gamma, ints)
    return ForwardDiff.derivative(x -> Dgamma2DR(X0, x, gamma, ints), beta)
end



function Riemann(X0, beta ,gamma, ints) 
    # beta = 1/T, gamma = -mu_B/T
    R_beta2 = Dbeta2DR(X0, beta, gamma, ints)
    R_gamma2 = Dgamma2DR(X0, beta, gamma, ints)
    R_beta_gamma = Dbeta_gammaDR(X0, beta, gamma, ints)
    R_beta2_gamma = Dbeta2_gammaDR(X0, beta, gamma, ints)
    R_gamma2_beta = Dgamma2_betaDR(X0, beta, gamma, ints)
    R_gamma3 = Dgamma3DR(X0, beta, gamma, ints)
    R_beta3 = Dbeta3DR(X0, beta, gamma, ints)

    mat_g = [R_beta2 R_beta_gamma; R_beta_gamma R_gamma2]
    g = det(mat_g)
    mat_R = [
        R_beta2 R_beta_gamma R_gamma2;
        R_beta3 R_beta2_gamma R_gamma2_beta;
        R_beta2_gamma R_gamma2_beta R_gamma3
    ]
    Riemann = 1/(2*g^2) * det(mat_R)
    return Riemann
end