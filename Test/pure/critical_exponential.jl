includet("../../src/pure/Rep.jl")
using DataFrames
using CSV
using BenchmarkTools
using Dierckx

function main1()
    # for alpha beta gamma

    Ts = range(131.015625, 130.01, length=20)
    firstline_path = "../../data/pure/1st.dat"
    df_in = CSV.read(firstline_path, DataFrame)
    df_unique = unique(df_in, :T)
    sort!(df_unique, :T)
    mu_unique = df_unique.muB
    T_unique = df_unique.T
    mu_of_T = Spline1D(T_unique, mu_unique, k=3, bc="error")

    ints = get_nodes(256;nodes2=256)
    data = zeros(length(Ts), 4)  # T, mu_B, C, chi_mumu
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        X0 = [-0.5, -0.5, -1.8, 0.80,0.80]  # 重置初始猜测
    
    
        mu_B = mu_of_T(T)   
        X0 = Tmu(T/hc, mu_B/hc, X0, ints)
        phi = X0[1:3]
        Phi1 = X0[4]
        Phi2 = X0[5]

        NewX = [phi..., Phi1, Phi2]
        C = DTTOmega(1.001*NewX, T/hc, mu_B/hc, ints) * T/hc 
        chi_mumu = Dmu2Omega(1.001*NewX, T/hc, mu_B/hc, ints)
        data[i, :] = [T, mu_B, C, chi_mumu]
  
    end
    df = DataFrame(data, [:T, :mu_B, :C, :chi_mumu])
    CSV.write("../../data/pure/CE1.dat", df)
end


function main2()
    # for alpha beta gamma

    Ts = range(131.01, 130.01, length=20)
    rhos = 0.01:0.01:3.00
    ints = get_nodes(256;nodes2=256)
    data = zeros(length(Ts), 4)  # T, mu_B, C, chi_mumu
    for (i, T) in enumerate(Ts)
        println("T = $T MEV")
        X0 = [-0.5, -0.5, -1.8, 0.80,0.80]  # 重置初始猜测
    
    
        mu_B = mu_of_T(T)   
        X0 = Tmu(T/hc, mu_B/hc, X0, ints)
        phi = X0[1:3]
        Phi1 = X0[4]
        Phi2 = X0[5]

        NewX = [phi..., Phi1, Phi2]
        C = DTTOmega(1.001*NewX, T/hc, mu_B/hc, ints) * T/hc 
        chi_mumu = Dmu2Omega(1.001*NewX, T/hc, mu_B/hc, ints)
        data[i, :] = [T, mu_B, C, chi_mumu]
  
    end
    df = DataFrame(data, [:T, :mu_B, :C, :chi_mumu])
    CSV.write("../../data/pure/CE2.dat", df)
end


function main3()
    # for alpha beta gamma

    T = 131.015625
    muBs = range(873.3582800185835+0.01, 873.3582800185835+0.1, length=30)

    ints = get_nodes(256;nodes2=256)
    data = zeros(length(muBs), 4)  # T, mu_B, C, chi_mumu
    for (i, muB) in enumerate(muBs)
        println("muB = $muB MEV")
        X0 = [-0.5, -0.5, -1.8, 0.80,0.80]  # 重置初始猜测
    
    

        X0 = Tmu(T/hc, muB/hc, X0, ints)
        phi = X0[1:3]
        Phi1 = X0[4]
        Phi2 = X0[5]

        NewX = [phi..., Phi1, Phi2]
        chi_mu = DmuOmega(1.001*NewX, T/hc, muB/hc, ints) / 0.16 #rho
        chi2mu = Dmu2Omega(1.001*NewX, T/hc, muB/hc, ints)
        data[i, :] = [T, muB, chi_mu, chi2mu]
  
    end
    df = DataFrame(data, [:T, :mu_B, :rho, :chi2mu])
    CSV.write("../../data/pure/CE3.dat", df)
end