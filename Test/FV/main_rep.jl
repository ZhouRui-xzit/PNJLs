include("../../src/Finite_Vol/Rep.jl")
using CSV 
using DataFrames

function main1()
    Ts = 300:-2:1
    mu_B = 800.0
   
    Rs = [3.0, 5.0, 10.0, 30.0]
    X0 = zeros(5)  # phi_u, phi_d, phi_s, Phi1, Phi2
    NewX = similar(X0)
    lens = length(Ts)
    data = zeros(lens, 5)  # T, mu_B, chi21, chi31, chi42
    for R in Rs 
        println("R = $R fm")
        ints = get_nodes(256, R)  # 每个R都要重新计算积分节点
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # 重置初始猜测值
        for (i, T) in enumerate(Ts)
            println("-------------------T = $T")
            NewX = Tmu(T/hc, mu_B/hc, X0, ints)
            flucts = Fluctuations(X0, T/hc, mu_B/hc, ints)
            data[i, :] = flucts
            X0 = NewX  # 使用上一个温度点的解作为下一个温度点的初始猜测值
            
        end
        df = DataFrame(data, [:T, :mu_B, :chi21, :chi31, :chi42])
        path = "../../data/FV/"
        CSV.write(path * "Fluctuations_mu=$(mu_B)_R=$(R).csv", df)
        println("Saved to $path/Fluctuations_mu=$(mu_B)_R=$(R).csv")
    end

end

main1()