using Revise
includet("../../src/magnetic/rep.jl")
using CSV 
using DataFrames


function main1()
    Ts = 300:-2:1
    mu_B = 800.0
    eBs = [0.20, 0.40, 0.60]
    X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # phi_u, phi_d, phi_s, Phi1, Phi2
    NewX = similar(X0)
    ints = get_nodes(128, 10)
    lens = length(Ts)
    data = zeros(lens, 5)  # T, mu_B, phi_u, phi_d, phi_s, Phi1, Phi2
    dir_path = "../../data/magnetic/"
    if !isdir(dir_path)
        println("Directory $dir_path does not exist, creating it.")
        mkpath(dir_path)
    end
    for eB in eBs
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # reset initial guess for each q
        for (i, T) in enumerate(Ts)
            println("T = $T")
            NewX = Tmu(T/hc, mu_B/hc, X0, eB*(1000/hc)^2, ints)
            data[i, :] =  numF(NewX, T/hc, mu_B/hc, eB*(1000/hc)^2, ints)

            X0 = NewX
        end
        df = DataFrame(data, [:T, :mu_B, :chi21, :chi31, :chi42])
        path = "../../data/magnetic/results_fluc_eB=$(eB).dat"
        CSV.write(path, df)
        println("Data saved to $path")
    end

end

main1()