include("TheroFunc.jl")
using CSV 
using DataFrames


function SolveOmega(X0, T, mu_B, ints)
    X0_typed = convert.(promote_type(eltype(X0), typeof(T), typeof(mu_B)), X0)
    
    fWrapper(Xs) = NewQuark_mu(Xs, T, mu_B, ints)

    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    NewX = res.zero
    phi = NewX[1:3]
    Phi1 = NewX[4]
    Phi2 = NewX[5]
    return -Omega(phi, Phi1, Phi2, T, mu_B, ints)
end

function DTOmega(X0, T, mu_B, ints)
    return ForwardDiff.derivative(x -> SolveOmega(X0, x, mu_B, ints), T)
end

function DTTOmega(X0, T, mu_B, ints)
    
    return ForwardDiff.derivative(x -> DTOmega(X0, x, mu_B, ints), T)
end

function DmuOmega(X0, T, mu_B, ints)
    return ForwardDiff.derivative(x -> SolveOmega(X0, T, x, ints), mu_B)
end

function DmumuOmega(X0, T, mu_B, ints)
    return ForwardDiff.derivative(x -> DmuOmega(X0, T, x, ints), mu_B)
end

function DmuTOmega(X0, T, mu_B, ints)
    return ForwardDiff.derivative(x -> DmuOmega(X0, x, mu_B, ints), T)
end




function Ther_Rep(NewX,T, mu_B, ints)

    chi_T = DTOmega(NewX, T, mu_B, ints)
    chi_TT = DTTOmega(NewX, T, mu_B, ints)
    chi_mu = DmuOmega(NewX, T, mu_B, ints)
    chi_mumu = DmumuOmega(NewX, T, mu_B, ints)
    chi_muT = DmuTOmega(NewX, T, mu_B, ints)

    P = SolveOmega(NewX, T, mu_B, ints)
    Cv = T * (chi_TT - chi_muT^2 / chi_mumu)
    
    Cp = T * (chi_TT - 2 * chi_T * chi_muT / chi_mu + (chi_T/chi_mu)^2 * chi_mumu)

    v_n_2 = (chi_T * chi_mumu - chi_mu * chi_muT) / (T * (chi_mumu * chi_TT - chi_muT^2))
    v_s_2 = (chi_T * chi_muT - chi_mu * chi_TT) / (mu_B * (chi_muT^2 - chi_TT * chi_mumu))
    E = T * chi_T + mu_B * chi_mu - P
    v_2 = (v_n_2 * T * chi_T + v_s_2 * mu_B * chi_mu) / (P + E)
    return [T*197.33, mu_B*197.33, chi_T, chi_TT, chi_mu, chi_mumu, chi_muT, P, Cv, Cp, v_2]
end

function zeta_and_eta()
    path = "Tmu0701.dat"
    my_data = CSV.read(path, DataFrame)
    df = DataFrame(my_data, [:T, :mu_B, :chi_T, :chi_TT, :chi_mu, :chi_mumu, :chi_Tmu, :P, :Cv, :Cp, :v_2])

    Ts = df.T 
    mu_Bs = df.mu_B
    S = df.chi_T
    rho = df.chi_mu
    chi_TT = df.chi_TT
    chi_mumu = df.chi_mumu
    chi_Tmu = df.chi_Tmu
    P = df.P
    Cv = df.Cv
    Cp = df.Cp
    cs = df.v_2
    E = Ts./hc .* S + mu_Bs./hc .* rho - P
    xi_alpha = 0.15656
    delta_Mu_xi = 20
    delta_T_xi = 20
    beta_delta = 1.5

    # f_new1 的实现
    term11 = ((abs.((mu_Bs ./ hc ./ 3 .- 873.42835 ./ 3 ./ hc) .* cos(xi_alpha) .- (Ts ./ hc .- 131.03 ./ hc) .* sin(xi_alpha))).^2) ./ (delta_Mu_xi ./ hc).^2
    term12 = ((abs.((mu_Bs ./ hc ./ 3 .- 873.42835 ./ 3 ./ hc) .* sin(xi_alpha) .+ (Ts ./ hc .- 131.03 ./ hc) .* cos(xi_alpha))).^(2 / beta_delta)) ./ (delta_T_xi ./ hc).^(2 / beta_delta)
    fnew1 = term11 .+ term12

    # f_new2 的实现（无除以 3 的版本）
    term21 = ((abs.((mu_Bs ./ hc .- 873.42835 ./ hc) .* cos(xi_alpha) .- (Ts ./ hc .- 131.03 ./ hc) .* sin(xi_alpha))).^2) ./ (delta_Mu_xi ./ hc).^2
    term22 = ((abs.((mu_Bs ./ hc .- 873.42835 ./ hc) .* sin(xi_alpha) .+ (Ts ./ hc .- 131.03 ./ hc) .* cos(xi_alpha))).^(2 / beta_delta)) ./ (delta_T_xi ./ hc).^(2 / beta_delta)
    fnew2 = term21 .+ term22

    # f_new3 的实现（PNJL 模型）
    term31 = ((abs.((mu_Bs ./ 3 ./ hc .- 873.42835 ./ 3 ./ hc) .* cos(xi_alpha) .- (Ts ./ hc .- 131.03 ./ hc) .* sin(xi_alpha))).^2) ./ (delta_Mu_xi ./ 3 ./ hc).^2
    term32 = ((abs.((mu_Bs ./ 3 ./ hc .- 873.42835 ./ 3 ./ hc) .* sin(xi_alpha) .+ (Ts ./ hc .- 131.03 ./ hc) .* cos(xi_alpha))).^(2 / beta_delta)) ./ (delta_T_xi ./ hc).^(2 / beta_delta)
    fnew3 = term31 .+ term32

    nu = 0.5

    # 计算 Xi1、Xi2 和 Xi3  关联长度
    Xi1 = (tanh.(fnew1) .* (1 .- 0.1^(2 / nu)) .+ 0.1^(2 / nu)).^(-nu / 2)
    Xi2 = (tanh.(fnew2) .* (1 .- 0.1.^(2/nu)) .+ 0.1.^(2/nu)).^(-nu/2)
    Xi3 = (tanh.(fnew3) .* (1 .- 0.1^(2 / nu)) .+ 0.1^(2 / nu)).^(-nu / 2)

    # 计算 eta_s 和 zeta_s 
    eta_s = Ts ./ hc .* Cp ./ (cs.^0.5) ./ (Xi2.^2) ./ S ./ Cv
    zeta_s = (cs.^0.5) .* Xi2 .* Cp .* (P .+ E) ./ S ./ Cv

    out = [Ts, mu_Bs, S, chi_TT, rho, chi_mumu, chi_Tmu, P, Cv, Cp, cs, Xi2,eta_s, zeta_s]
    df_out = DataFrame(out, [:T, :mu_B, :S, :chi_TT, :rho, :chi_mumu, :chi_Tmu, :P, :Cv, :Cp, :cs, :Xi2,:eta_s,:zeta_s])
    out_path = "Tmu_more0701.dat"
    CSV.write(out_path, df_out, writeheader=true)
end