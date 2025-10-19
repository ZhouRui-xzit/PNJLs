############################################################
# PNJL Omega(T, mu_B) 计算核心（含 minimal 方案）
#
# 约定（重要）：
# 1) 单位：内部统一使用 fm 与 fm^-1。
#    - 传入 T, mu_B 必须是 fm^-1（若外部以 MeV 给定，请在调用处除以 hc）。
# 2) 本文件实现 minimal 方案：
#    - 真空项：固定区间 [0, Lambda_f] 上的固定 GL 网格积分。
#    - 热  项：在 [0, pmax] 上做“过渡层窗口”分段，窗口段用较密节点，其余段用稀节点；
#              节点为固定 GL 网格，且基础节点只在全局构造一次，积分时仅线性映射。
#    - 不做全域自适应，避免频繁重建网格，保证 Trho 100ms 量级性能，Tmu 10ms 量级性能。
# 3) 保留 integrate_GL_adapt 以便对照或临时回退；默认 Omega 使用 minimal 实现。


############################################################


############################################################
# ForwardDiff 使用注意事项（本文件已按此修正）：
# - 任何由 phi/T/mu_B 推导出的中间容器，元素类型必须随输入标量
#   漂移（可能是 Dual）。禁止使用 Vector{Float64}、Tuple{Float64,...} 等
#   固定浮点容器；一律用类型参数化 + promote_type。
# - 常数尽量用 zero(S)/one(S)/S(数值) 构造。
# - 仅在确知“不需要梯度”的路径上，才显式 value() 剥离 Dual。
############################################################

include("constant.jl")

using FastGaussQuadrature
using ForwardDiff
using NLsolve
using NaNMath
# ========================= minimal 方案参数 =========================
# 固定节点数（窗口内/外）
const N_WIN       = 64          # 窗口段较密节点
const N_OUT       = 32          # 非窗口段较稀节点
# 过渡层窗口宽度系数与最小半宽（fm^-1）
const K_WIN       = 8.0         # Delta = K_WIN * T * pF / EF
const DELTA_MIN   = 0.03        # 防止 T 很小时窗口退化为零
# 端点合并阈值（避免产生碎片段）
const EPS_MERGE   = 1e-3        # fm^-1
# 热项 p 上限（fm^-1）
const PMAX_TH     = 20.0

# 预先缓存两套基础 GL 节点（在 (-1,1) 上）
const GL_X_OUT, GL_W_OUT = gausslegendre(N_OUT)
const GL_X_WIN, GL_W_WIN = gausslegendre(N_WIN)

# =============== 可选：自适应 GL（保留以便对照或回退） ===============
function integrate_GL_adapt(g::Function; a::Real=0.0, b::Real=1.0,
                            atol::Real=1e-12, rtol::Real=1e-10,
                            N0::Int=128, Nmax::Int=512)
    @assert N0 >= 2 "N0 must be >= 2"
    @assert Nmax >= N0 "Nmax must be >= N0"

    gl(N) = begin
        x, w = gausslegendre(N)          # (-1,1)
        p  = (b - a)/2 .* x .+ (b + a)/2 # 映射到 [a,b]
        wp = (b - a)/2 .* w
        sum(wp .* g.(p))
    end

    Iprev = gl(N0)
    N = N0
    while true
        N2 = min(2N, Nmax)
        I  = gl(N2)
        if abs(I - Iprev) <= max(atol, rtol * max(abs(I), abs(Iprev), 1.0)) || N2 == Nmax
            return I
        end
        N, Iprev = N2, I
    end
end

# =========================== 物理子模块 ===========================
function AA(x, T, Phi1, Phi2)
    # A = 1 + 3*Phi1*e^{-x/T} + 3*Phi2*e^{-2x/T} + e^{-3x/T}
    term1 = exp(-x / T)
    term2 = exp(-2.0 * x / T)
    term3 = exp(-3.0 * x / T)
    return 1.0 + 3.0 * Phi1 * term1 + 3.0 * Phi2 * term2 + term3
end

function AAbar(x, T, Phi1, Phi2)
    # Abar = 1 + 3*Phi2*e^{-x/T} + 3*Phi1*e^{-2x/T} + e^{-3x/T}
    term1 = exp(-x / T)
    term2 = exp(-2.0 * x / T)
    term3 = exp(-3.0 * x / T)
    return 1.0 + 3.0 * Phi2 * term1 + 3.0 * Phi1 * term2 + term3
end



function chiral(phi, theta)
    @inbounds begin
        sigma = phi[1:3]
        eta   = phi[4:6]
        sigma_u, sigma_d, sigma_s = sigma
        eta_u,   eta_d,   eta_s   = eta
    end

    # 二次项
    term1 = 2 * G_f * (sum(sigma .^ 2) + sum(eta .^ 2))

    # 't Hooft 六次项（含 theta 耦合）
    term2 = 4 * K_f * (
          cos(theta) * (sigma_u * sigma_d * sigma_s)
        + sin(theta) * (eta_u   * eta_d   * eta_s)
        - cos(theta) * (sigma_u * eta_d   * eta_s + eta_u * eta_d * sigma_s + eta_u * eta_s * sigma_d)
        - sin(theta) * (sigma_u * sigma_s * eta_d + sigma_u * sigma_d * eta_s + sigma_d * sigma_s * eta_u)
    )

    return term1 + term2
end




function Mass(phi, theta)
    """计算三种夸克的有效质量"""
    sigma = phi[1:3]
    eta   = phi[4:6]

    mu_s = m0[1] + 4 * G_f * sigma[1] + 2 * K_f * (
            cos(theta) * (sigma[2] * sigma[3] - eta[2] * eta[3]) -
            sin(theta) * (sigma[2] * eta[3] + eta[2] * sigma[3])
    )

    mu_ps = 4 * G_f * eta[1] - 2 * K_f * (
            cos(theta) * (sigma[2] * eta[3] + eta[2] * sigma[3]) +
            sin(theta) * (sigma[2] * sigma[3] - eta[2] * eta[3])
    )
    
    mass_u = sqrt(mu_s^2 + mu_ps^2)

    md_s = m0[2] + 4 * G_f * sigma[2] + 2 * K_f * (
            cos(theta) * (sigma[3] * sigma[1] - eta[3] * eta[1]) -
            sin(theta) * (sigma[3] * eta[1] + eta[3] * sigma[1])
    )
    
    md_ps = 4 * G_f * eta[2] - 2 * K_f * (
            cos(theta) * (sigma[3] * eta[1] + eta[3] * sigma[1]) +
            sin(theta) * (sigma[3] * sigma[1] - eta[3] * eta[1])
    )
    
    mass_d = sqrt(md_s^2 + md_ps^2)
    
    ms_s = m0[3] + 4 * G_f * sigma[3] + 2 * K_f * (
            cos(theta) * (sigma[1] * sigma[2] - eta[1] * eta[2]) -
            sin(theta) * (sigma[1] * eta[2] + eta[1] * sigma[2])
    )
    
    ms_ps = 4 * G_f * eta[3] - 2 * K_f * (
            cos(theta) * (sigma[1] * eta[2] + eta[1] * sigma[2]) +
            sin(theta) * (sigma[1] * sigma[2] - eta[1] * eta[2])
    )
    
    mass_s = sqrt(ms_s^2 + ms_ps^2)
    

    return [mass_u, mass_d, mass_s]
end

E_of_p(p, mass) = sqrt(p*p + mass*mass)

function calc_U(T, Phi1, Phi2)
    # 对数型 Polyakov 势
    Ta = a0 + a1 * (T0/T) + a2 * (T0/T)^2
    Tb = b3 * (T0 / T)^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    U = T^4 * (-0.5 * Ta * Phi2 * Phi1 + Tb * log(value))
    return U
end

# =========================== integrand ===========================
@inline function vacand(p, phi, theta)
    # 真空 integrand（包含相空间因子 p^2/(2*pi^2)）
    M  = Mass(phi, theta)
    sE = zero(eltype(M))
    @inbounds for f in 1:3
        sE += sqrt(p*p + M[f]*M[f])
    end
    return (-6 * sE) * (p*p) / (2*pi^2)
end

@inline function thand(p, phi, Phi1, Phi2, T, muB, theta)
    # 热 integrand（包含相空间因子 p^2/(2*pi^2)）
    M   = Mass(phi, theta)
    mu3 = muB/3
    slog = zero(eltype(M))
    @inbounds for f in 1:3
        Ef = sqrt(p*p + M[f]*M[f])
        slog += log(AA(Ef - mu3, T, Phi1, Phi2)) + log(AAbar(Ef + mu3, T, Phi1, Phi2))
    end
    return (-2*T*slog) * (p*p) / (2*pi^2)
end

# ===================== minimal: 分段与固定网格 =====================
# 合并端点：把间距 < tol 的相邻端点合并，避免碎片段
function _merge_knots!(knots::AbstractVector{S}; tol::Real=EPS_MERGE) where {S<:Real}
    sort!(knots)
    lens = length(knots)
    j = 1
    @inbounds for i in 2:lens
        if knots[i] - knots[j] > S(tol)
            j += 1
            knots[j] = knots[i]
        end
    end
    resize!(knots, j)
    return knots
end

# 构造热项分段：按每个味的 Fermi 动量在 [0, pmax] 上切出窗口段
function _make_th_segments(phi, T, mu_B, theta; pmax::Real=PMAX_TH,
                           k::Real=K_WIN, delta_min::Real=DELTA_MIN,
                           tol::Real=EPS_MERGE)
    # 统一标量类型：可能是 Float64 或 Dual
    S = promote_type(eltype(phi), typeof(T), typeof(mu_B))

    M   = Mass(phi, theta)                        # eltype(M) 与 eltype(phi) 一致
    mu3 = S(mu_B)/S(3)

    windows = NTuple{2,S}[]                # 窗口端点容器（类型随 S）
    @inbounds for f in 1:3
        Mf = S(M[f])
        # 相对阈值：接近等于时不进窗口（防止 μ≈M 处误入）
        tolμ = S(1e-12) * (abs(mu3) + abs(Mf) + one(S))
    if mu3 > Mf + tolμ
        # 正部截断，避免 sqrt(负数)：若 μ3≈Mf 但差值被浮点/AD 扰动成负，就设为 0
        Δ2 = mu3*mu3 - Mf*Mf
        Δ2 = ifelse(Δ2 > zero(S), Δ2, zero(S))
        if Δ2 == zero(S)
            # 差值极小：相当于无费米面，直接跳过窗口
            # （也可选择构造一个极窄窗口，但通常跳过更稳）
        else
            pF = sqrt(Δ2)
            EF = sqrt(pF*pF + Mf*Mf)             # 或 hypot(pF, Mf)（若你确认对 Dual OK）
            # Δ = k * T * pF / EF，EF>0 已由分支保证；为稳妥再做一次最小值夹逼
            den = ifelse(EF > zero(S), EF, one(S))
            delta = S(k) * S(T) * pF / den
            delta = ifelse(delta < S(delta_min), S(delta_min), delta)

            a = max(zero(S), pF - delta)
            b = min(S(pmax), pF + delta)
            if b > a + S(1e-12)
                push!(windows, (a, b))
            end
        end
    end
    end

    # 端点集：0, pmax, 以及各窗口端点
    knots = S[zero(S), S(pmax)]
    for (a,b) in windows
        push!(knots, a); push!(knots, b)
    end
    _merge_knots!(knots; tol=tol)

    # 分段并标记是否在窗口内
    segs = Vector{Tuple{S,S,Bool}}()
    @inbounds for i in 1:length(knots)-1
        a = knots[i]; b = knots[i+1]
        if b <= a + S(1e-15)
            continue
        end
        mid = (a + b) / S(2)
        iswin = any(mid >= aw && mid <= bw for (aw,bw) in windows)
        push!(segs, (a, b, iswin))
    end
    return segs
end

# 在 [a,b] 上用缓存基础节点做一次固定 GL 积分（不含相空间外的权改动）
@inline function _gl_sum_segment(a, b, iswin::Bool, f::Function)
    # 让节点映射后的 p 与权 wp 也随标量类型 S 漂移
    S = promote_type(typeof(a), typeof(b))
    if iswin
        x, w = GL_X_WIN, GL_W_WIN
    else
        x, w = GL_X_OUT, GL_W_OUT
    end
    half = (S(b) - S(a)) / S(2)
    mid  = (S(b) + S(a)) / S(2)
    p  = half .* S.(x) .+ mid
    wp = half .* S.(w)
    return sum(wp .* f.(p))
end

# -------------------- 真空项（单段） --------------------
@inline function _integrate_vac_minimal(phi, theta; Λf::Real=Lambda_f)
    S = promote_type(eltype(phi), Float64)
    a, b = zero(S), S(Λf)
    x, w = GL_X_OUT, GL_W_OUT
    half = (b - a) / S(2)
    mid  = (b + a) / S(2)
    p  = half .* S.(x) .+ mid
    wp = half .* S.(w)
    return sum(wp .* vacand.(p, Re(phi), theta))
end

# -------------------- 热项（窗口分段 + 固定网格） --------------------
function _integrate_th_minimal(phi, Phi1, Phi2, T, mu_B, theta;
                               pmax::Real=PMAX_TH,
                               k::Real=K_WIN, delta_min::Real=DELTA_MIN,
                               tol::Real=EPS_MERGE)
    segs  = _make_th_segments(phi, T, mu_B, theta; pmax=pmax, k=k, delta_min=delta_min, tol=tol)
    total = zero(promote_type(eltype(phi), typeof(T), typeof(mu_B)))
    @inbounds for (a,b,iswin) in segs
        total += _gl_sum_segment(a, b, iswin, p -> thand(p, phi, Phi1, Phi2, T, mu_B, theta))
    end
    return total
end
# ============================== Omega ==============================
function Omega(phi, Phi1, Phi2, T, mu_B, theta;
               Λf::Real=Lambda_f, pmax::Real=PMAX_TH,
               # 自适应参数仍保留供回退时使用，但 minimal 分支不使用：
               atol::Real=1e-12, rtol::Real=1e-10, N0::Int=128, Nmax::Int=512)
    # 物理局部项
    chi = chiral(phi, theta);
    U   = calc_U(T, Phi1, Phi2)

    # minimal 方案：分段 + 固定缓存网格
    Iv  = _integrate_vac_minimal(phi, theta; Λf=Λf)
    Ith = _integrate_th_minimal(phi, Phi1, Phi2, T, mu_B, theta; pmax=pmax)
    
    return chi + U + Iv + Ith
end

# 如需临时回退到“自适应版本”的 Omega，可改为：
# function Omega(...)
#     chi = chiral(phi); U = calc_U(T, Phi1, Phi2)
#     Iv  = integrate_GL_adapt(p -> vacand(p, phi); a=0.0, b=Λf, atol=atol, rtol=rtol, N0=N0, Nmax=Nmax)
#     Ith = integrate_GL_adapt(p -> thand(p, phi, Phi1, Phi2, T, mu_B); a=0.0, b=pmax, atol=atol, rtol=rtol, N0=N0, Nmax=Nmax)
#     return chi + U + Iv + Ith
# end

# =========================== dOmega 导数 ===========================
function dOmega_dphi(phi, Phi1, Phi2, T, mu_B, theta)
    return ForwardDiff.gradient(x -> Omega(x, Phi1, Phi2, T, mu_B, theta), phi)
end

function dOmega_dPhi1(phi, Phi1, Phi2, T, mu_B, theta)
    return ForwardDiff.derivative(x -> Omega(phi, x, Phi2, T, mu_B, theta), Phi1)
end

function dOmega_dPhi2(phi, Phi1, Phi2, T, mu_B, theta)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, x, T, mu_B, theta), Phi2)
end

function dOmega_dT(phi, Phi1, Phi2, T, mu_B, theta)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, x, mu_B, theta), T)
end

function dOmega_dmu_B(phi, Phi1, Phi2, T, mu_B, theta)
    return ForwardDiff.derivative(x -> Omega(phi, Phi1, Phi2, T, x, theta), mu_B)
end

# =========================== 方程与求解 ===========================
function function_Quark_core(x_array, T, mu_B, theta)
    # 核心 8 维驻点条件
    phi = x_array[1:6]
    Phi1 = x_array[7]
    Phi2 = x_array[8]
    dphi  = dOmega_dphi(phi, Phi1, Phi2, T, mu_B, theta)
    dP1   = dOmega_dPhi1(phi, Phi1, Phi2, T, mu_B, theta)
    dP2   = dOmega_dPhi2(phi, Phi1, Phi2, T, mu_B, theta)
    return [dphi[1], dphi[2], dphi[3], dphi[4], dphi[5], dphi[6], dP1, dP2]
end

function Quark_mu(x_array, T, mu_B, theta)
    # T-mu_B 模式：8 维非线性方程
    T_out = promote_type(eltype(x_array), typeof(T), typeof(mu_B), typeof(theta))
    fvec = zeros(T_out, 8)
    fvec .= function_Quark_core(x_array, T, mu_B, theta)
    return fvec
end

function Quark_rho(x_array, T, rho_B, theta)
    # T-rho_B 模式：11 维非线性方程（含密度约束）
    T_out = promote_type(eltype(x_array), typeof(T), typeof(rho_B), typeof(theta))
    fvec = zeros(T_out, 11)

    X0   = x_array[1:8]
    phi  = x_array[1:6]
    Phi1 = x_array[7]
    Phi2 = x_array[8]

    mu_B = x_array[9] * 3
    fvec[1:8] .= function_Quark_core(X0, T, mu_B, theta)

    rho_now = -1 .* dOmega_dmu_B(phi, Phi1, Phi2, T, mu_B, theta)
    fvec[9] = x_array[9] - x_array[10]
    fvec[10] = x_array[10] - x_array[11]
    fvec[11] = rho_now / rho0 - rho_B
    return fvec
end

function Tmu(T, mu_B, X0, theta)
    # 给定 T, mu_B（均为 fm^-1）求解 8 维驻点
    fWrapper(Xs) = Quark_mu(Xs, T, mu_B, theta)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    return res.zero
end

function Trho(T, rho_B, X0, theta)
    # 给定 T, rho_B 求解 11 维驻点
    fWrapper(Xs) = Quark_rho(Xs, T, rho_B, theta)
    res = nlsolve(fWrapper, X0, autodiff=:forward)
    return res.zero
end
