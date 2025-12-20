using Revise
includet("some_main_el.jl")
includet("../../src/Finite_Vol/QuarkPhase.jl")
includet("derv.jl")


function quark_rhos()
    R = 20.0
    es = [0.0, 0.3]
    CEP_data = CSV.read("../../data/FV/CEP_result_el.csv", DataFrame)
    for (R, e) in Iterators.product([R], es)
        T_CEP = CEP_data[(CEP_data.R .== R) .& (CEP_data.e .== e), :T_CEP][1]
        println("R=$R e=$e T_CEP=$T_CEP")
        main_Trho(R, e, T_CEP)
        QP_Trho(R, e)
    end
end


function quark_mus()
    R = 20.0
    es = [0.0, 0.3]
    
    for (R, e) in Iterators.product([R], es)
        CEP_data = CSV.read("../../data/FV/Trho_Maxwell_R=$(R)_e=$(e).csv", DataFrame)
        T_CEP = CEP_data[end, :T]
        muc = CEP_data[end, :mu_star]
        
        println("R=$R, e=$e: T_CEP=$T_CEP MeV, μc=$muc MeV")
        main_Tmu(R, e, muc)
        QP_Tmu(R, e, T_CEP)
    end
end

function quark_mus_without_CEP()
    R = 20.0
    es = [0.7, 1.0]
    
    for (R, e) in Iterators.product([R], es)
        #CEP_data = CSV.read("../../data/FV/Trho_Maxwell_R=$(R)_e=$(e).csv", DataFrame)
        #T_CEP = CEP_data[end, :T]
        #muc = CEP_data[end, :mu_star]
        T_CEP = 0.0
        muc = 900.0
        println("R=$R, e=$e: T_CEP=$T_CEP MeV, μc=$muc MeV")
        main_Tmu(R, e, muc)
        QP_Tmu(R, e, T_CEP)
    end
end



function test_quark_mus()
    R = 20.0
    es = [0.7]
    
    for (R, e) in Iterators.product([R], es)
        #CEP_data = CSV.read("../../data/FV/Trho_Maxwell_R=$(R)_e=$(e).csv", DataFrame)
        #T_CEP = CEP_data[end, :T]
        #muc = CEP_data[end, :mu_star]
        T_min = 10.0
        T_max = 40.0
        
        println("R=$R, e=$e: T_min=$T_min MeV, T_max=$T_max MeV")
        main_muT(R, e, T_min, T_max)
        find_muT(R, e)
    end
end

function find_muT(R, e)
    df = CSV.read("../../data/FV/mu_B_T_R=$(R)_e=$(e).csv", DataFrame)
    
    Ts_unique = unique(df.T)
    sort!(Ts_unique)
    
    mu_maxs = zeros(Float64, length(Ts_unique))
    dphi_maxs = zeros(Float64, length(Ts_unique))
    
    for (i, T) in enumerate(Ts_unique)
        # 筛选当前温度的数据
        df_T = df[df.T .== T, :]
        sort!(df_T, :mu_B)
        
        mu_Bs = Float64.(df_T.mu_B)
        phius = Float64.(df_T.phi_u)
        
        # 计算导数
        dphius = zeros(Float64, length(phius))
        derivation!(phius, mu_Bs, dphius)
        

        # 找极大值
        inds_max, vals_max = findmaxima(dphius)
    if !isempty(inds_max)
        # 选导数最大的极大值对应的位置
        k = inds_max[argmax(vals_max)]
        mu_peak = mu_Bs[k]
        dphi_peak = dphius[k]
        println("global max at mu_B=$mu_peak, value=$dphi_peak")
        mu_maxs[i] = mu_peak
        dphi_maxs[i] = dphi_peak
        else
            mu_maxs[i] = NaN
            dphi_maxs[i] = NaN
            @warn "未找到极大值: T=$T MeV"
        end
    end
    
    # 保存结果
    df_out = DataFrame(T=Ts_unique, mu_B_max=mu_maxs, dphi_dmu_max=dphi_maxs)
    outpath = "../../data/FV/muT_maxima_R=$(R)_e=$(e).csv"
    CSV.write(outpath, df_out)
    println("结果已保存至 $outpath")
    
    return df_out
end