using Revise

includet("../../src/Finite_Vol/pnjl_FV.jl")
using DataFrames
using CSV
using BenchmarkTools
using Dates


function main_Tmu(paras)
    println("Time:", Dates.now())
    Ts = 220:-1:100.0
    muc = 296.40677563816973 * 3
    mu_B1 = range(0.0,muc, length=30)
    #mu_B2 = 301.0*3:0.5:304.14*3
    mu_B = vcat(mu_B1)
    a, b, c = paras
    ints = get_nodes_el(128, a, b, c)
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    lens = length(Ts) * length(mu_B)
    data = zeros(lens, 7)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, MU) in enumerate(mu_B)
        println("mu_B = $MU MEV")
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # 重置初始猜测
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            
            X0 = Tmu(T/hc, MU/hc, X0, ints)
            data[idx, :] = [T, MU, X0...]
        end
    end
    df = DataFrame(data, [:T, :mu_B, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write("../../data/FV/T_mu_B_scan_el=$a.dat", df)
end

function main_Trho(paras)
    println("Time:", Dates.now())
    #T_CEP =  129.960   #a=b = 30 , c = 30
    #T_CEP =  129.007   #a=b = 10 , c = 30
    #T_CEP =  128.29   #a=b = 7 , c = 30
    #T_CEP =  127.28  #a=b = 5 , c = 30
    #T_CEP = 124.68 #a=b = 3 , c = 30
    #T_CEP = 121.0  #a=b = 2 , c = 30
    T_CEP = 106.08  #a=b = 1 , c = 30
    

    T1s = T_CEP:-0.01:T_CEP-0.1
    T2s = (T_CEP-0.1)-0.1:-0.1:T_CEP-1.0
    T3s = T_CEP-2.0:-2.0:10.0
    Ts = vcat(T1s, T2s, T3s)

    a, b, c = paras

    rho1s = 3.00:-0.01:0.01
    rho2s = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    rho2s = reverse(rho2s)
    rhos = vcat(rho1s, rho2s)
    
    X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 320/hc, 320/hc, 320/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    ints1 = get_nodes_el(128, a, b, c)
    ints2 = get_nodes_el(256, a, b, c)
    ints3 = get_nodes_el(512, a, b, c)

    lens = length(Ts) * length(rhos)
    data = zeros(lens, 9)  # T, rho_B, mu_B/3, P, phi_u, phi_d, phi_s, Phi1, Phi2
   
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        if T<= 50.0
            ints = ints2
        elseif T<= 20.0
            ints = ints3
        else
            ints = ints1
        end

        X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 320/hc, 320/hc, 320/hc]  # 重置初始猜测
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
    CSV.write("../../data/FV/T_rho_B_scan_el=$a.dat", df)
end



function main_Tmu_el()
    Rs = [30.0, 10.0, 7.0, 5.0, 3.0, 2.0, 1.5, 0.9]
    mu_B = 0.0
    Ts = 300.0:-2.0:10.0

    lens = length(Ts) * length(Rs)
    data = zeros(lens, 7)  # T, R, phi_u, phi_d, phi_s, Phi1, Phi2

    for (i, R) in enumerate(Rs)
        println("R = $R fm")
        ratio = 4.0
        theta = pi / ratio  # 圆形横向截面
        #a = abs(R * cos(theta))
        #b = abs(R * sin(theta))
        a=R; b=R;
        c = 30.0
        ints = get_nodes_el(128, a, b, c)
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # 重置初始猜测
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            
            X0 = Tmu(T/hc, mu_B/hc, X0, ints)
            data[idx, :] = [T, R, X0...]
        end
    end

    outpath = "../../data/FV/T_mu0_el_N.csv"
    df = DataFrame(data, [:T, :R, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write(outpath, df)
    println("结果已保存至 $outpath")
end


function main_Tmu_el()
 
    mu_B = 0.0
    Ts = 300.0:-2.0:10.0
    as = [1.0, 3.0, 5.0, 7.0, 10.0, 30.0]
    lens = length(Ts) * length(as)
    data = zeros(lens, 7)  # T, R, phi_u, phi_d, phi_s, Phi1, Phi2
    
    b = 1.0
    c = 30.0
    for (i, a) in enumerate(as)
        println("a = $a fm")
        ints = get_nodes_el(128, a, b, c)
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # 重置初始猜测
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            
            X0 = Tmu(T/hc, mu_B/hc, X0, ints)
            data[idx, :] = [T, a, X0...]
        end
    end

    outpath = "../../data/FV/T_mu0_el_vara_N.csv"
    df = DataFrame(data, [:T, :a, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
    CSV.write(outpath, df)
    println("结果已保存至 $outpath")
end








if abspath(PROGRAM_FILE) == @__FILE__
    paras = [10.0, 10.0,]
    main_Trho(paras)
end