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


using StaticArrays
# 基本几何量（按给定参数化）
# r_θ, r_φ, r_θθ, r_θφ, r_φφ
rθ(a,b,c, θ,φ)  = SVector(a*cos(θ)*cos(φ),  b*cos(θ)*sin(φ), -c*sin(θ))
rφ(a,b,c, θ,φ)  = SVector(-a*sin(θ)*sin(φ), b*sin(θ)*cos(φ),  0.0)
rθθ(a,b,c,θ,φ)  = SVector(-a*sin(θ)*cos(φ), -b*sin(θ)*sin(φ), -c*cos(θ))
rθφ(a,b,c,θ,φ)  = SVector(-a*cos(θ)*sin(φ),  b*cos(θ)*cos(φ),  0.0)
rφφ(a,b,c,θ,φ)  = SVector(-a*sin(θ)*cos(φ), -b*sin(θ)*sin(φ),  0.0)

# 第一基本形式系数
function first_fundamental(a,b,c, θ,φ)
    u = rθ(a,b,c,θ,φ); 
    v = rφ(a,b,c,θ,φ)
    E = dot(u,u)
    F = dot(u,v)
    G = dot(v,v)
    return E,F,G,u,v
end

# 单位外法向
function normal_and_area_factor(a,b,c, θ,φ)
    u = rθ(a,b,c,θ,φ); 
    v = rφ(a,b,c,θ,φ)
    w = cross(u,v)
    J = norm(w)                      # |r_θ × r_φ| = √(EG - F^2)
    n̂ = w / J                       # 取外法向
    return n̂, J
end

# 第二基本形式系数
function second_fundamental(a,b,c, θ,φ, n̂)
    e = dot(rθθ(a,b,c,θ,φ), n̂)
    f = dot(rθφ(a,b,c,θ,φ), n̂)
    g = dot(rφφ(a,b,c,θ,φ), n̂)
    return e,f,g
end

# 标量场：面素、H、K
function local_invariants(a,b,c, θ,φ)
    E,F,G,_,_ = first_fundamental(a,b,c,θ,φ)
    n̂, J     = normal_and_area_factor(a,b,c,θ,φ)
    e,f,g     = second_fundamental(a,b,c,θ,φ,n̂)
    denom     = E*G - F^2
    H = (E*g - 2F*f + G*e) / (2*denom)     # H = (κ1+κ2)/2 约定
    K = (e*g - f^2) / denom
    return J, H, K
end

function integrate_over_ellipsoid(a, b, c; rtol=1e-9)
    # φ 积分，返回三个分量的积分
    function φint(θ)
        quadgk(0.0, 2π; rtol=rtol) do φ
            J, H, K = local_invariants(a, b, c, θ, φ)
            SVector(J, H*J, K*J)   # 用SVector替代元组
        end |> x -> x[1]           # x[1] 是 SVector(∫J, ∫HJ, ∫KJ)
    end

    # θ 积分，分别对三个分量积分
    result = quadgk(0.0, π; rtol=rtol) do θ
        φint(θ)
    end |> x -> x[1]               # x[1] 是 SVector(∫∫J, ∫∫HJ, ∫∫KJ)

    return result[1], result[2], result[3]
end
