using Revise

includet("../../src/Finite_Vol/pnjl_FV.jl")
includet("derv.jl")
using DataFrames
using CSV
using BenchmarkTools
using Dates
using Peaks


function main_Tmu(R, e, muc)
    println("Time:", Dates.now())

    Ts = 220:-1:10.0
   
    mu_B1 = range(0.0, muc, length=30)
    #mu_B2 = 850.0:2.0:900.0
    mu_B = unique(vcat(mu_B1))
    a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)

    ints = get_nodes_el(128, a, b, c, modes="D")
    X00 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    lens = length(Ts) * length(mu_B)
    data = zeros(lens, 7)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, MU) in enumerate(mu_B)
        println("mu_B = $MU MEV")
        X0 = X00
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            
            X0 = Tmu(T/hc, MU/hc, X0, ints)
            data[idx, :] = [T, MU, X0...]
            if j==1
                X00 = X0
            end
        end
    end
    df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/FV/T_mu_B_R=$(R)_e=$(e).csv", df)
end

function main_Trho(R, e, T_CEP)
 
    a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)

    modes = "D"
    #T_CEP = 107.921875 # R=30.0 e=0.0
    #T_CEP = 98.9765625 # R=30.0 e=0.3
    #T_CEP = 81.2734375 # R=30.0 e=0.7
    #T_CEP = 69.5546875 # R=30.0 e=1.0



    T1s = T_CEP:-0.02:T_CEP-0.1
    T2s = (T_CEP-0.1)-0.1:-0.2:T_CEP-1.0
    T3s = range(T_CEP-1.0, T_CEP-10.0, length=5)  
    T4s = range(T_CEP-10.0, 10.0, length=30)       
    Ts = unique(vcat(T1s, T2s, T3s, T4s))


    rho1s = 3.00:-0.01:0.01
    rho2s = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
   
    rhos = sort(vcat(rho1s, rho2s), rev=true)
    
    X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 1200/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    ints = get_nodes_el(300, a, b, c;modes=modes)
    

    lens = length(Ts) * length(rhos)
    data = zeros(lens, 9)  # T, rho_B, mu_B/3, P, phi_u, phi_d, phi_s, Phi1, Phi2
   
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        

        X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 1200/hc]  # 重置初始猜测
        for (j, rho_B) in enumerate(rhos)
            idx = (i-1)*length(rhos) + j
            
            X0 = Trho(T/hc, rho_B, X0, ints)
            mu = X0[6] * hc 
            phi = X0[1:3] 
            Phi1 = X0[4]
            Phi2 = X0[5]
            P = -Omega(phi, Phi1, Phi2, T/hc, mu/hc, ints) 
            data[idx, :] = [T, rho_B, mu, P, X0[1:5]...]
        end
    end
    df = DataFrame(data, [:T, :rho_B, :mu, :P, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/FV/T_rho_B_R=$(R)_e=$(e).csv", df)
end



function main_Tmu_es()
    println("Time:", Dates.now())

    Ts = 300.0:-0.1:10.0
    mu_B = 0.0
    R = 20.0
    es = 0.0:0.1:1.0
    phius = similar(Ts)
    phi1s = similar(Ts)
    dphius = similar(Ts)
    dphi1s = similar(Ts)
    max_phius = similar(es)
    max_phi1s = similar(es)
    for (i, e) in enumerate(es)
        println("e = $e ")
        a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)
        ints = get_nodes_el(128, a, b, c, modes="D")
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
        for (j, T) in enumerate(Ts)
            #println("e = $e , T = $T MEV")
            X0 = Tmu(T/hc, mu_B/hc, X0, ints)
            phius[j] = X0[1]
            phi1s[j] = X0[4]
        end
        derivation!(phius, Ts, dphius)
        derivation!(phi1s, Ts, dphi1s)
        inds_max1, _ = findmaxima(dphius)
        inds_max2, _ = findmaxima(dphi1s)
        max_phius[i] = Ts[inds_max1[1]]
        max_phi1s[i] = Ts[inds_max2[1]]
    end
    df = DataFrame(e=es, T_phi_u_max=max_phius, T_Phi1_max=max_phi1s)
    CSV.write("T_mu_B_R=$(R)_el_maxima.dat", df)
end



function main_Tmu_Rs()
    println("Time:", Dates.now())

    Ts = 300.0:-0.1:10.0
    mu_B = 0.0
    Rs = 100.0:-5.0:20.0
    e = 0.0
    phius = similar(Ts)
    phi1s = similar(Ts)
    dphius = similar(Ts)
    dphi1s = similar(Ts)
    max_phius = similar(Rs)
    max_phi1s = similar(Rs)
    for (i, R) in enumerate(Rs)
        println("R = $R fm ")
        a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)
        ints = get_nodes_el(128, a, b, c, modes="D")
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
        for (j, T) in enumerate(Ts)
            #println("e = $e , T = $T MEV")
            X0 = Tmu(T/hc, mu_B/hc, X0, ints)
            phius[j] = X0[1]
            phi1s[j] = X0[4]
        end
        derivation!(phius, Ts, dphius)
        derivation!(phi1s, Ts, dphi1s)
        inds_max1, _ = findmaxima(dphius)
        inds_max2, _ = findmaxima(dphi1s)
        max_phius[i] = Ts[inds_max1[1]]
        max_phi1s[i] = Ts[inds_max2[1]]
    end
    df = DataFrame(R=Rs, T_phi_u_max=max_phius, T_Phi1_max=max_phi1s)
    CSV.write("T_mu_B_e=$(e)_maxima.dat", df)
end


function parametrize_deformation(R, δ;para=2.0,scale=1.0)
    """
    δ: 变形参数 (0 ≤ δ < ∞)
    para: 调节变形幅度的参数
    scale: +1 压扁 -1 拉长
    - δ = 0: 球形 (a=b=c=R)
    - δ > 0 (且 scale=1.0): 扁平椭球(a=b > c)
    - δ > 0 (且 scale=-1.0): 拉长椭球(a=b < c)
    - 表面积单调递增
    """
    V = (4/3)*π*R^3
    
    # 基于 β₂ 的简化
    β₂ = tanh(δ)  # 保证 β₂ < 1

    a = R * (1 + para * β₂)^(scale * 2/3)
    b = a
    c = V / ((4/3)*π*a*b)
    
    return a, b, c
end


function main_muT(R, e, T_CEP)
    println("Time:", Dates.now())

    mu_Bs = 0.0:10.0:1200.0
    Ts = range(200.0, T_CEP+1, length=30)


    a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)

    ints = get_nodes_el(128, a, b, c, modes="D")
    X00 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    lens = length(Ts) * length(mu_Bs)
    data = zeros(lens, 7)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        X0 = X00
        for (j, mu_B) in enumerate(mu_Bs)
            idx = (i-1)*length(mu_Bs) + j
            
            X0 = Tmu(T/hc, mu_B/hc, X0, ints)
            data[idx, :] = [T, mu_B, X0...]
            if j==1 
                X00 = X0
            end
        end
    end
    df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/FV/T_mu_B_R=$(R)_e=$(e).csv", df)
end