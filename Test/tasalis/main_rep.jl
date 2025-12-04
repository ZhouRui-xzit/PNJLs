include("../../src/Tsallis/rep.jl")
using Plots
using LaTeXStrings
using DataFrames, CSV

function main_rep(;mu_B=10.0)
   
    Ts = 300.0:-2.0:10.0
    #mu_B = 10.0
    eBs = [0.13, 0.30]
    qs = [1.0000001, 1.01, 1.05]
    Nodes = get_nodes(128, 20)
    nodes = gauleg(0, 1, 128)
    sol = zeros(length(Ts), 6)
    data = zeros(length(Ts), 5)
    P0 = 20.1235
    sol[:, 1] = Ts
    for eB in eBs
        println("Calculations for eB = ", eB, ) 
        for q in qs
            println("q = ", q)
            X0 = [-0.01, -0.01, -0.20, 0.8, 0.8]
            for (i, T) in enumerate(Ts)
                if i % 10 == 0
                    println("---------------------------Calculating T = ", T, " MeV; Step ", i, " of ", length(Ts))
                end
                NewX = Tmu(T/hc, mu_B/hc, X0, eB*(1000/hc)^2, q, Nodes, nodes)
                
                sol[i, 2:6] = NewX
                data[i, :] = Ther_Rep(NewX*1.0001, T/hc, mu_B/hc, eB*(1000/hc)^2, q, Nodes, nodes, P0)
                X0 = NewX
            end
            df1 = DataFrame(sol, [:T, :phi_u, :phi_d, :sigma, :Phi1, :Phi2])
            df2 = DataFrame(data, [:T,:mu_B,:P,:E,:cs2])
            CSV.write("sol_eB=$(eB)_q=$(q).csv", df1)
            CSV.write("Rep_eB=$(eB)_q=$(q).csv", df2)
        end
    end 
    
end



