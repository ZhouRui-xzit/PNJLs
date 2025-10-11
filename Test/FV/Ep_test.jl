using FastGaussQuadrature

function total_mean_curvature_ellipsoid(a,b,c; nθ=200, nφ=400)
    x, wθ = gausslegendre(nθ)                 # x ∈ [-1,1] = cosθ
    φ  = range(0, 2π, length=nφ+1)[1:end-1]
    wφ = 2π/nφ
    acc = 0.0
    for (i, xi) in pairs(x)
        s2 = 1 - xi^2                         # sin^2θ
        row = 0.0
        for ϕ in φ
            row += sqrt(a^2*s2*cos(ϕ)^2 + b^2*s2*sin(ϕ)^2 + c^2*xi^2)
        end
        acc += row * wθ[i] * wφ               # 注意：这里没有 sinθ
    end
    return 2 * acc                             # 不是 4！
end

# 轴对称（一维积分更快更稳）
function total_mean_curvature_spheroid(a,c; n=400)  # a=a, b=a, c=c
    # θ∈[0,π/2] Gauss-Legendre
    x, w = gausslegendre(n)             # x∈[-1,1]
    θ = (x .+ 1) .* (π/4)               # 映射到 [0,π/2]
    jac = π/4
    sθ = sin.(θ); cθ = cos.(θ)
    f  = sqrt.(a^2 .* sθ.^2 .+ c^2 .* cθ.^2) .* sθ
    return 8π * sum(f .* w) * jac
end
