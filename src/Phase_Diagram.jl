include("constants.jl")
include("TheroFunc.jl")
include("Functions_Quark.jl")
include("Rep.jl")

using NLsolve
using ForwardDiff

function Tmu(T, mu_B, X0, ints)
    Xquark = X0
    T = T / hc
    mu_B = mu_B / hc

    
    fWrapper(Xs) = NewQuark_mu(Xs, T, mu_B, ints)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    NewX = res.zero
    
    return NewX
end

function Trho(T, rho_B, X0, ints)
    
    T = T / hc

    fWrapper(Xs) = Quark_rho(Xs, T, rho_B, ints)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    NewX = res.zero
    return NewX
end


