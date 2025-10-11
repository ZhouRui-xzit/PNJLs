include("../../src/pure/Rep.jl")
using CSV 
using DataFrames

function main1()
    Ts = 300:-2:50
    mu_B = 800.0
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    NewX = similar(X0)
    lens = length(Ts)
    data = zeros(lens, 5)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, T) in enumerate(Ts)
        println("T = $T")
        NewX = Tmu(T/hc, mu_B/hc, X0)
        data[i, :] =  Fluctuations(X0, T/hc, mu_B/hc)
        X0 = NewX
    end
    df = DataFrame(data, [:T, :mu_B, :chi21, :chi31, :chi42])
    CSV.write("../../data/pure/Fluctuations_mu=$mu_B.csv", df)
    println("Saved to ../../data/pure/Fluctuations_mu=$mu_B.csv")
end

main1()