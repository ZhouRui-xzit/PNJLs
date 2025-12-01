using Revise
includet("../../src/magnetic/pnjl_magnetic.jl")
using Plots
using LaTeXStrings
using DataFrames, CSV


function main_eBs()
    Ts = 305:-2.0:1.0
    X0 = zeros(5)
    NewX = zeros(5)
    mu_B = 10.0
    eBs = [0.20, 0.40, 0.60 ,0.80]
    Nodes1 = get_nodes(128, 100)
    nodes = gauleg(0.0, 100, 500)
    lens = length(Ts) * length(eBs)
    data = zeros(lens, 7)
    
    for (i, eB) in enumerate(eBs)
        X0 = [-0.01, -0.01, -0.20, 0.8, 0.8]
        println("eB = ", eB)
        for (j, T) in enumerate(Ts)
            println("-------------------------------T = ", T)
            idx = (i - 1) * length(Ts) + j
            NewX = Tmu(T/hc, mu_B/hc, X0, eB*(1000/hc)^2, Nodes, nodes)
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


function generate_rhos2(n::Int)
    rhos2 = Float64[]
    for e in -n:-2
        push!(rhos2, 1.0 * 10.0^e)
        push!(rhos2, 5.0 * 10.0^e)
    end
    #push!(rhos2, 0.05)  # 确保包含终点
    return sort(rhos2)
end


function main_rho(;T=50.0)
    rhos1 = 4.00:-0.01:0.01
    rhos2 = generate_rhos2(13)
    rhos2 = reverse(rhos2)
    rhos = vcat(rhos1, rhos2)
    rhos = sort(rhos, rev=true)
    #rhos = 1.00:-0.01:0.01
    eBs = [0.13]
    mub = 1200.11978  # MeV
    X0 = [-0.15025, -0.14023, -2.06331, 0.05312, 0.07734, mub/hc]
   
    data = zeros(length(rhos), 9) 
    Nodes = get_nodes(500, 20)
    nodes = gauleg(0, 1, 256)
    for eB in eBs
        println("eB = ", eB, " T = ", T)
        for (i, rho) in enumerate(rhos)
            println("-------------------------------rho = ", rho)
            NewX = Trho(T/hc, rho, X0, eB*(1000/hc)^2, Nodes, nodes)
            data[i, 1] = T  # MeV
            data[i, 2] = rho  # rho_B/rho
            data[i, 3] = eB # GeV^2
            data[i, 4:8] = NewX[1:5] # phi_u, phi_d, phi_s, Phi1, Phi2
            data[i, 9] = NewX[6] * hc  # mu_B MeV
            X0 = NewX 
            println("mu_B = ", data[i, 9])
        end
        outpath = "../../data/magnetic/Trho_eB_$(round(eB,digits=2))_T=$(T).dat"
        df = DataFrame(data, [:T, :rho_B, :eB, :phiu, :phid, :phis, :Phi1, :Phi2, :mu_B])
        CSV.write(outpath, df)
        println("Data saved to ", outpath)
    end
end


#@time main_rho(30)