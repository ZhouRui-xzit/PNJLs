include("constants.jl")

using FastGaussQuadrature 


function gauleg(a, b, n)
    x, w = gausslegendre(n)
    x_mapped = (b - a) / 2 .* x .+ (b + a) / 2
    w_mapped = (b - a) / 2 .* w
    return x_mapped, w_mapped
end

function get_nodes(p_num)
    p1, w1 = gauleg(0.0,Lambda_f, p_num)
    p2, w2 = gauleg(0.0,20.0, p_num)
    int1 = [p1,w1]
    int2 = [p2,w2]
    return int1, int2
end

function  int_same(int)
    p, w = int
    int_same = w .* p.^2 ./ (2*pi^2)
    return int_same
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