using Revise

includet("../../src/Finite_Vol/pnjl_FV.jl")
using DataFrames
using CSV
using BenchmarkTools
using Dates


function main_Tmu(;R=30.0,e=0.0)
    println("Time:", Dates.now())

    Ts = 220:-1:75
    muc = 299.409884100486 * 3
    mu_B1 = range(0.0, muc, length=30)

    mu_B = vcat(mu_B1)
    a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)

    ints = get_nodes_el(128, a, b, c, modes="D")
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
    CSV.write("../../data/FV/T_mu_B_scan_el=$e.dat", df)
end

function main_Trho(;R=30.0,e=0.0)
 
    a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)

    modes = "D"
    T_CEP =  107.921875

    T1s = T_CEP:-0.02:T_CEP-0.1
    T2s = (T_CEP-0.1)-0.1:-0.2:T_CEP-1.0
    T3s = range(T_CEP-1.0, T_CEP-10.0, length=5)  # 修改：使用 range
    T4s = range(T_CEP-10.0, 10.0, length=30)       # 修改：使用 range
    Ts = unique(vcat(T1s, T2s, T3s, T4s))
    Ts = 107:-1:10
    T_test = [20.0]

    rho1s = 3.00:-0.01:0.01
    rho2s = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
   
    rhos = sort(vcat(rho1s, rho2s))
    
    X0 = [-1.8,-1.8, -2.2, 0.01,0.01, 320/hc, 320/hc, 320/hc]  # phi_u, phi_d, phi_s, Phi1, Phi2, muB_div3, aux1, aux2
    ints1 = get_nodes_el(200, a, b, c;modes=modes)
    ints2 = get_nodes_el(400, a, b, c;modes=modes)
    ints3 = get_nodes_el(800, a, b, c;modes=modes)

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
    CSV.write("../../data/FV/T_rho_B_single_el=$e.dat", df)
end



function main_Tmu_el()
    Rs = [30.0, 10.0, 5.0, 1.5]
    mu_B = 0.0
    Ts = 300.0:-2.0:10.0

    lens = length(Ts) * length(Rs)
    data = zeros(lens, 7)  # T, R, phi_u, phi_d, phi_s, Phi1, Phi2

    for (i, R) in enumerate(Rs)
        println("R = $R fm")
        ratio = 4.0
        theta = pi / ratio  # 圆形横向截面

        a=R; b=R;
        c = 30.0
        ints = get_nodes_el(128, a, b, c; modes="N")
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
  #@  para = 1.8
    a = R * (1 + para * β₂)^(scale * 2/3)
    b = a
    c = V / ((4/3)*π*a*b)
    
    return a, b, c
end




function main_Tmu_equal_V()
    mu_B = 0.0
    Ts = 300:-2.0:10.0
    R = 30.0
    #V = (4/3)*pi*R^3
    delta_s = [0.0, 0.3, 0.7, 1.0]
    lens = length(Ts) * length(delta_s)
    data = zeros(lens, 7)  # T, R, phi_u, phi_d, phi_s, Phi1, Phi2
    for (i, e) in enumerate(delta_s)
        println("e = $e ")
        a, b, c = parametrize_deformation(R, e;para=3.0,scale=-1.0)
        ints = get_nodes_el(128, a, b, c, modes="D")
        X0 = [-0.01, -0.01, -0.40, 0.8, 0.8]  # 重置初始猜测
        for (j, T) in enumerate(Ts)
            idx = (i-1)*length(Ts) + j
            X0 = Tmu(T/hc, mu_B/hc, X0, ints)
            data[idx, :] = [T, e, X0...]

        end

    end
outpath = "../../data/FV/D_V_R=$(R).csv"
df = DataFrame(data, [:T, :e, :phi_u, :phi_d, :phi_s, :Phi1, :Phi2])
CSV.write(outpath, df)
println("结果已保存至 $outpath")
    

end

