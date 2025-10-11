include("../../src/magnetic/pnjl_magnetic.jl")
using Plots
using LaTeXStrings
using DataFrames, CSV


function main_eBs()
    Ts = 305:-2.0:1.0
    X0 = zeros(5)
    NewX = zeros(5)
    mu_B = 10.0
    eBs = [0.20, 0.40, 0.60 ,0.80]
    Nodes = get_nodes(128, 100)

    lens = length(Ts) * length(eBs)
    data = zeros(lens, 7)
    
    for (i, eB) in enumerate(eBs)
        X0 = [-0.01, -0.01, -0.20, 0.8, 0.8]
        println("eB = ", eB)
        for (j, T) in enumerate(Ts)
            println("-------------------------------T = ", T)
            idx = (i - 1) * length(Ts) + j
            NewX = Tmu(T/hc, mu_B/hc, X0, eB*(1000/hc)^2, Nodes)
            data[idx, 1] = T
            data[idx, 2] = eB
            data[idx, 3:end] = NewX
            X0 = NewX
        end
    end
    outpath = "../../data/magnetic/Tmu_all_eB.dat"
    df = DataFrame(data, [:T, :eB, :phi1, :phi2, :phi3, :Phi1, :Phi2])
    CSV.write(outpath, df)
    println("Data saved to ", outpath)
end

function main_eBs2()
    Ts = 305:-2.0:1.0
    X0 = zeros(5)
    NewX = zeros(5)
    mu_B = 0.0
    eBs = [0.20, 0.40, 0.60 ,0.80]
    Nodes = get_nodes(128, 10)

    lens = length(Ts) * length(eBs)
    data = zeros(lens, 5)
    phi0 = [1.84321051058777, 1.84321051058777, 2.226907587459626]
    mpi = 135 / hc  # pi 介子质量
    fpi = 87.9 / hc # pi 介子衰变常数
    masses = [5.5, 5.5, 140.7] ./ hc  # 当前夸克质量
    for (i, eB) in enumerate(eBs)
        X0 = [-0.01, -0.01, -0.20, 0.8, 0.8]
        println("eB = ", eB)
        for (j, T) in enumerate(Ts)
            println("-------------------------------T = ", T)
            idx = (i - 1) * length(Ts) + j
            NewX = Tmu(T/hc, mu_B/hc, X0, eB*(1000/hc)^2, Nodes)
            data[idx, 1] = T
            data[idx, 2] = eB
            phi = NewX[1:3]
            
            sigmau = (2*masses[1] )/(mpi^2 * fpi^2) * (-phi[1]-phi0[1]) + 1
            sigmad = (2*masses[2] )/(mpi^2 * fpi^2) * (-phi[2]-phi0[2]) + 1
            sigmas = (2*masses[3] )/(mpi^2 * fpi^2) * (-phi[3]-phi0[3]) + 1
            data[idx, 3] = sigmau
            data[idx, 4] = sigmad
            data[idx, 5] = sigmas
            X0 = NewX
        end
    end
    outpath = "../../data/magnetic/Tmu_all_eB.dat"
    df = DataFrame(data, [:T, :eB, :sigmau, :sigmad, :sigmas])
    CSV.write(outpath, df)
    println("Data saved to ", outpath)
end



function main_rho(T)
    #rhos1 = 5:-0.01:0.01
    #rhos2 = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
    #rhos2 = reverse(rhos2)
    #rhos = vcat(rhos1, rhos2)
    #rhos = sort(rhos)
    rhos = 0.01:0.01:5.0
    eBs = [0.05]
    X0 = [-2.8,-2.8, -2.2, 0.01,0.01, 320/hc, 320/hc, 320/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    data = zeros(length(rhos), 8) 
    Nodes = get_nodes(128, 100)
    for eB in eBs
        println("eB = ", eB, " T = ", T)
        for (i, rho) in enumerate(rhos)
            println("-------------------------------rho = ", rho)
            NewX = Trho(T/hc, rho, X0, eB*(1000/hc)^2, Nodes)
            data[i, 1] = rho  # rho_B/rho
            data[i, 2] = eB # GeV^2
            data[i, 3:7] = NewX[1:5] # phi_u, phi_d, phi_s, Phi1, Phi2
            data[i, 8] = NewX[6] * 3 * hc # mu_B in MeV
            X0 = NewX 
            println("mu_B = ", data[i, 8])
        end
        outpath = "../../data/magnetic/Trho_eB_$(round(eB,digits=2)).dat"
        df = DataFrame(data, [:rho_B, :mu_B, :eB])
        CSV.write(outpath, df)
        println("Data saved to ", outpath)
    end
end


main_eBs2()

#@time main_rho(30)