include("constants.jl")
include("gauleg.jl")

include("TheroFunc.jl")
include("Functions_Quark.jl")


using CSV, DataFrames
using DelimitedFiles #读写数据





function Trho(T, rho_B, X0, ints)
    
    T = T / hc

    NewX, success = Root_rho(X0, T, rho_B, ints)
    if success == false
        println("Root finding failed")
        data = [T / hc, 0.0, 0.0, 0.0]
        return X0,data
    end
    phi = NewX[1:3]
    Phi1 = NewX[4]
    Phi2 = NewX[5]
    mu_u = NewX[6]
    mu_B = mu_u * 3
    P = -Omega(phi, Phi1, Phi2, T,mu_B, ints)
    rho_B_norm = -dOmgea_dmu_B / rho0

    data = [T*hc, mu_u*hc, rho_B_norm, P]
    return NewX, data
end


function main(T_start, T_end)
    inte1, inte2 = get_nodes(128)
    int_same1 = int_same(inte1)
    int_same2 = int_same(inte2)
    ints = [inte1, inte2, int_same1, int_same2]

    T_loop = T_end:-1:T_start
    rho_loop = 3.00:-0.01:0.10

    for T = T_loop
        X0 = [-1.8, -1.8, -2.1, 0.8, 0.8, 320 / hc, 320 / hc, 320 / hc]
        for rho_B = rho_loop
            X0, data = Trho(T, rho_B, X0, ints)
        end
    end
    
end