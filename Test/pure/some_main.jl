include("../../src/pure/Pnjl_pure.jl")
using DataFrames
using CSV
using BenchmarkTools



function main_Tmu()
    Ts = 200:-0.01:100.0
    mu_B1 = range(0.0, 100.0*3, length=10)
    mu_B2 = range(100.0*3, 280.0*3, length=20)
    mu_B3 = range(280.0*3, 290.0*3, length=5)

    mu_Bs = unique!(vcat(collect(mu_B1), collect(mu_B2), collect(mu_B3)))
    mu_Bs = [291.0*3]
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
    CSV.write("../../data/pure/T_mu_B_scan.dat", df)
end



function main_Trho()

    ints1 = get_nodes(200, nodes2=200)
    ints2 = get_nodes(256, nodes2=500)

    T1s = range(131.015625, 131.01, length=10)
    T2s = range(131.00, 130.00, length=10)
    T3s = range(130.0, 10.0, length=50)
    Ts = unique!(vcat(collect(T1s), collect(T2s), collect(T3s)))

    rho1s = 3.00:-0.01:0.01
    rho2s = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    rho2s = reverse(rho2s)
    rhos = vcat(rho1s, rho2s)
    X00 = [-1.8,-1.8, -2.2, 0.01,0.01, 1000/hc]  # 重置初始猜测
    X0 = similar(X00)
    lens = length(Ts) * length(rhos)
    data = zeros(lens, 8)  # T, rho_B, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        if T >= 30.0
            ints = ints1
        else
            ints = ints2
        end
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
    rhos1 = 3.00:-0.01:0.01
    rhos2 = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    rhos2 = reverse(rhos2)
    rhos = vcat(rhos1, rhos2)
    ints = get_nodes(128, nodes2=128)
    X0 = [-0.9,-0.9, -1.5, 0.01,0.01, 320/hc, 320/hc, 320/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    data = zeros(length(rhos), 3) 
    for (i, rho_B) in enumerate(rhos)
        println("T = $T, rho_B = $rho_B")
        X0 = Trho(T/hc, rho_B, X0, ints)
        #mu = X0[6] * hc 
        data[i, :] = [T, rho_B, mu]
    end
    df = DataFrame(data, [:T, :rho_B, :mu_B])
    CSV.write("../../data/pure/T=$T.dat", df)
end

function single_mu()
    T = 300.0
    mu_B = 0.0
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    X_sol = Tmu(T/hc, mu_B/hc, X0)
    #println("T = $T, mu_B = $mu_B")
    #println("Solution: phi_u=$(X_sol[1]), phi_d=$(X_sol[2]), phi_s=$(X_sol[3]), Phi1=$(X_sol[4]), Phi2=$(X_sol[5])")
end

function single_rho()
    T = 10.0
    rho_B = 3.0
    X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 320/hc, 320/hc, 320/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    X_sol = Trho(T/hc, rho_B, X0)
    mu = X_sol[6] * hc 
    #println("T = $T, rho_B = $rho_B")
    #println("Solution: phi_u=$(X_sol[1]), phi_d=$(X_sol[2]), phi_s=$(X_sol[3]), Phi1=$(X_sol[4]), Phi2=$(X_sol[5]), mu_B=$mu")
end
