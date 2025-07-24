include("Tools.jl")
include("Phase_Diagram.jl")
include("constants.jl")

using DelimitedFiles  # 引入标准库，无需额外安装
using DataFrames
using CSV

function main()
    i = 1
    for mu = mu_loop
        X0 = [-1.8, -1.8, -2.2, 0.0038, 0.0038]
        #X0 = [-0.5, -0.5, -1.8, 0.8, 0.8]
        for T = T_loop
            NewX = Tmu(T, mu, X0, ints)
            XSS[i] = NewX
            i = i + 1
            println("T=", T)
        end
    end

    j = 1
    for mu = mu_loop
        for T = T_loop
            X0 = XSS[j]
            my_data[j, :] = Ther_Rep(X0, T/hc, mu/hc, ints)
            j = j + 1
            println("T=", T*hc)
        end
    end

        # 将数据写入文件
    df = DataFrame(my_data, [:T, :mu_B, :chi_T, :chi_TT, :chi_mu, :chi_mumu, :chi_Tmu, :P, :Cv, :Cp, :v_2])
    CSV.write(output_file, df, writeheader=true)
    println("Data saved to ", output_file)

end



int1, int2 = get_nodes(128)
int_same1 = int_same(int1)
int_same2 = int_same(int2)
ints = [int1, int2, int_same1, int_same2]


T_loop = 210.0:-1:130.0


mu_loop = 0.0:1:5.0
lens = length(T_loop) * length(mu_loop)

my_data = zeros(lens, 11)
XSS = Vector{Vector{Float64}}(undef, lens)
# 定义输出文件路径
output_file = "Tmu0701.dat"


