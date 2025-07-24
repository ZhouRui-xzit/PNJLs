
include("constants.jl")
include("Tools.jl")
include("TheroFunc.jl")

using NLsolve
using ForwardDiff




function NewQuark_mu(x_array, T, mu_B, int_params)
    # 动态创建与输入类型兼容的数组
    T_out = promote_type(eltype(x_array), typeof(T), typeof(mu_B))
    fvec = zeros(T_out, 5)

    fvec[1], fvec[2], fvec[3], fvec[4], fvec[5] = function_Quark_core(x_array, T, mu_B, int_params)
    return fvec
end




function Quark_rho(x_array, T, rho_B, int_params)
    T_out = promote_type(eltype(x_array), typeof(T), typeof(rho_B))

    fvec = zeros(T_out, 8)
    X0 = x_array[1:5]
    phi = x_array[1:3]
    Phi1 = x_array[4]
    Phi2 = x_array[5]

    mu_B = x_array[6] * 3
    fvec[1], fvec[2], fvec[3], fvec[4], fvec[5] = function_Quark_core(X0, T, mu_B, int_params)
    rho_now = -1 .* dOmgea_dmu_B(phi, Phi1, Phi2, T, mu_B, int_params) 

    fvec[6] = x_array[6] - x_array[7]
    fvec[7] = x_array[7] - x_array[8]
    fvec[8] = rho_now / rho0 - rho_B
    return fvec
end


