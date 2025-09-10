include("constants.jl")

using FastGaussQuadrature 
using MeshGrid

function gauleg(a, b, n)
    x, w = gausslegendre(n)
    x_mapped = (b - a) / 2 .* x .+ (b + a) / 2
    w_mapped = (b - a) / 2 .* w
    return x_mapped, w_mapped
end

function get_nodes(p_num, n_num)
    p1, w1 = gauleg(-20.0, 20.0, p_num)
    ns = range(0, n_num-1)
    wn = ones(length(ns))


    P1, N1 = meshgrid(p1, ns)
    W1, WN = meshgrid(w1, wn)
    Alpha = ones(n_num)
    for i in 2:n_num
        Alpha[i] = 2
    end
    W = W1 .* WN .* Alpha ./ (4 * pi^2)
    ints = [P1, N1, W]
    return ints
end



function dAdB(A, B)
    n = length(A)
    dAdB = similar(A)
    
    # 中心差分（内部点）
    for i in 2:n-1
        dAdB[i] = (B[i+1] - B[i-1]) / (A[i+1] - A[i-1])
    end
    
    # 边界处理（前向/后向差分）
    dAdB[1] = (B[2] - B[1]) / (A[2] - A[1])
    dAdB[end] = (B[end] - B[end-1]) / (A[end] - A[end-1])
    
    return dAdB
end