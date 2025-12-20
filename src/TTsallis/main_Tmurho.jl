using Revise

includet("../../src/TTsallis/Pnjl_TT.jl")

using DataFrames
using CSV
using BenchmarkTools
using Dates
using Peaks


function main_Tmu()
    println("Time:", Dates.now())
    q = 1.001
    Ts = 220:-0.1:100.0
    muc = 893.42785
    mu_B1 = range(0.0, 150*3, length=10)
    mu_B2 = range(150*3, muc, length=20)
    mu_B = unique(vcat(mu_B1, mu_B2))

    ints = get_nodes(128, nodes2=128)
    X00 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    lens = length(Ts) * length(mu_B)
    data = zeros(lens, 8)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, MU) in enumerate(mu_B)
        println("mu_B = $MU MEV")
        X0 = X00
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            
            X0 = Tmu(T/hc, MU/hc, X0, q, ints)
            rho_B = -1 .* dOmega_dmu_B(X0[1:3], X0[4], X0[5], T/hc, MU/hc, q, ints) / rho0
            data[idx, :] = [T, MU, rho_B, X0...]
            if j==1
                X00 = X0
            end
        end
    end
    df = DataFrame(data, [:T, :mu_B, :rho_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("T_mu_B_scan.csv", df)
end

function main_Trho(q)
 
    TCEP_ab =  112.87451171875
    ints = get_nodes(300, nodes2=300)
    
    T1s = range(TCEP_ab, 112.80, length=5)
    T2s = range(112.80, 110.0, length=10)
    T3s = range(110.0, 10.0, length=20)



    Ts = unique(vcat(T1s, T2s, T3s))

    Ts = [TCEP_ab]
    rho1s = 3.00:-0.005:0.01
    rho2s = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
   
    rhos = sort(vcat(rho1s, rho2s), rev=true)
    
    X00 = [-1.8,-1.8, -2.2, 0.01,0.01, 1200/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB
    lens = length(Ts) * length(rhos)
    data = zeros(lens, 9)  # T, rho_B, mu_B P, phi_u, phi_d, phi_s, Phi1, Phi2
   
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        

        X0 = X00
        for (j, rho_B) in enumerate(rhos)
            idx = (i-1)*length(rhos) + j
            
            X0 = Trho(T/hc, rho_B, X0, q, ints)
            mu = X0[6] * hc 
            phi = X0[1:3] 
            Phi1 = X0[4]
            Phi2 = X0[5]
            P = -Omega(phi, Phi1, Phi2, T/hc, mu/hc, q, ints) 
            data[idx, :] = [T, rho_B, mu, P, X0[1:5]...]
            if j==1
                X00 = X0
            end
        end
    end
    df = DataFrame(data, [:T, :rho_B, :mu, :P, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("T_rho_above.csv", df)
end






