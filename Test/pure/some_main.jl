include("../../src/pure/Pnjl_pure.jl")
using DataFrames
using CSV
using BenchmarkTools



function main_Tmu()
    Ts = 220:-1.0:120.0
    mu_B1 = range(0.0, 100.0*3, length=10)
    mu_B2 = range(100.0*3, 873.50, length=30)


    mu_Bs = unique!(vcat(collect(mu_B1), collect(mu_B2)))

    ints = get_nodes(128;nodes2=128)
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    lens = length(Ts) * length(mu_Bs)
    data = zeros(lens, 8)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, MU) in enumerate(mu_Bs)
        println("mu_B = $MU MEV")
        #if 
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # 重置初始猜测
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            
            X0 = Tmu(T/hc, MU/hc, X0, ints)
            NewX = 1.0001 .* X0
            rho = -dOmega_dmu_B(NewX[1:3], NewX[4], NewX[5], T/hc, MU/hc, ints) / 0.16 
            data[idx, :] = [T, MU, X0...,rho]
        end
    end
    df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2, :rho_B])
    CSV.write("../../data/pure/T_mu_B_scan.csv", df)
end



function main_Trho()

   

    T1s = range(131.03375244140625, 131.03, length=10)
    T2s = range(131.00, 130.00, length=10)
    T3s = range(130.0, 10.0, length=30)
    Ts = unique!(vcat(collect(T1s), collect(T2s), collect(T3s)))

    rho1s = 3.00:-0.01:0.01
    rho2s = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    rho2s = reverse(rho2s)
    rhos = vcat(rho1s, rho2s)

    ints = get_nodes(300;nodes2=300)
    X00 = [-1.8, -1.8, -2.2, 0.01,0.01, 1000/hc]  # 重置初始猜测
    X0 = similar(X00)
    lens = length(Ts) * length(rhos)
    data = zeros(lens, 8)  # T, rho_B, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        X0 = X00
        for (j, rho_B) in enumerate(rhos)
            idx = (i-1)*length(rhos) + j
            
            X0 = Trho(T/hc, rho_B, X0, ints)
            mu = X0[6] * hc 
            
            data[idx, :] = [T, rho_B, mu, X0[1:5]...]
            if j==1
                X00 = X0
            end
        end
    end
    df = DataFrame(data, [:T, :rho_B, :mu, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/pure/T_rho_B_scan1123.dat", df)
end





function main_mu(mu_B)
    Ts = 300:-1:0.1
    
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    data = zeros(length(Ts), 7)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, T) in enumerate(Ts)
        println("T = $T, mu_B = $mu_B")
        X0 = Tmu(T/hc, mu_B/hc, X0)
        data[i, :] = [T, mu_B, X0...]
    end
    df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/pure/mu_B=$mu_B.dat", df)
end

function main_rho(T)
    rhos = 3.00:-0.01:0.50

    ints = get_nodes(128, nodes2=128)
    X0 = [-0.9,-0.9, -1.5, 0.01,0.01, 1000/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    data = zeros(length(rhos), 3) 
    for (i, rho_B) in enumerate(rhos)
        println("T = $T, rho_B = $rho_B")
        X0 = Trho(T/hc, rho_B, X0, ints)
        mu = X0[6] * hc 
        data[i, :] = [T, rho_B, mu]
    end
    df = DataFrame(data, [:T, :rho_B, :mu_B])
    CSV.write("../../data/pure/T=$T.dat", df)
end

