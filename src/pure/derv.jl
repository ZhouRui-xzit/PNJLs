"""
HCHIR: given N data points (X[i],Y[i]), computes a minimax (Remez-like) polynomial
approximation of degree (M-1) in power basis.

Returns A of length (M+1):
- A[1:M] are coefficients for p(x)=A1 + A2 x + ... + A[M] x^(M-1)
- A[M+1] stores the (nonnegative) max error amplitude (or negative on early return,
  matching the Fortran behavior)
"""
function hchir(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}, N::Int, M::Int)
    @assert length(X) == N && length(Y) == N
    M  = min(M, N-1, 19)
    M1 = M + 1

    A  = zeros(Float64, M1)
    IX = ones(Int, M1)
    H  = zeros(Float64, M1)

    HA = 0.0

    IX[1]  = 1
    IX[M1] = N
    L = (N - 1) ÷ M
    J = L
    for i in 2:M
        IX[i] = J + 1
        J += L
    end

    while true
        # label 20 block in Fortran
        HH = 1.0
        @inbounds for i in 1:M1
            A[i] = Float64(Y[IX[i]])
            H[i] = -HH
            HH   = -HH
        end

        # label 50 block
        @inbounds for j in 1:M
            II = M1
            Y2 = A[II]
            H2 = H[II]
            for i in j:M
                D  = Float64(X[IX[II]]) - Float64(X[IX[M1 - i]])
                Y1 = A[M - i + j]
                H1 = H[M - i + j]
                A[II] = (Y2 - Y1) / D
                H[II] = (H2 - H1) / D
                II = M - i + j
                Y2 = Y1
                H2 = H1
            end
        end

        HH = -A[M1] / H[M1]
        @inbounds for i in 1:M1
            A[i] = A[i] + H[i] * HH
        end

        # label 80 block: convert to power basis
        @inbounds for j in 1:(M-1)
            II = M - j
            D  = Float64(X[IX[II]])
            Y2 = A[II]
            for k in (M1 - j):M
                Y1 = A[k]
                A[II] = Y2 - D * Y1
                Y2 = Y1
                II = k
            end
        end

        HM = abs(HH)
        if HM <= HA
            A[M1] = -HM
            return A
        end
        A[M1] = HM
        HA = HM

        IM = IX[1]
        H1 = HH
        j = 1

        # label 100 block: find max error point not in IX
        @inbounds for i in 1:N
            if j <= M1 && i == IX[j]
                if j < M1
                    j += 1
                end
            else
                # Horner evaluate p(x)
                H2 = A[M]
                for k in (M-1):-1:1
                    H2 = H2 * Float64(X[i]) + A[k]
                end
                H2 -= Float64(Y[i])
                if abs(H2) > HM
                    HM = abs(H2)
                    H1 = H2
                    IM = i
                end
            end
        end

        if IM == IX[1]
            return A
        end

        # locate insertion position I
        I = 1
        while I <= M1 && IM >= IX[I]
            I += 1
        end
        if I > M1
            I = M1
        end

        H2 = (isodd(I) ? -HH : HH)

        if H1 * H2 >= 0.0
            IX[I] = IM
            continue
        end

        if IM < IX[1]
            for jj in M:-1:1
                IX[jj+1] = IX[jj]
            end
            IX[1] = IM
            continue
        end

        if IM > IX[M1]
            for jj in 2:M1
                IX[jj-1] = IX[jj]
            end
            IX[M1] = IM
            continue
        end

        IX[I-1] = IM
    end
end


"""
derivation!(a,b,c): array derivative algorithm (Fortran derivation clone)

- a, b of length n, with b being the grid (not necessarily uniform)
- writes derivative into c (Float64), same length

Interior points i=2..n-1: uses the given 3-point nonuniform formula.
Endpoints: uses HCHIR on 5-point windows to fit cubic to c and extrapolate to b[1], b[n].
"""
function derivation!(a::AbstractVector{<:Real},
                     b::AbstractVector{<:Real},
                     c::AbstractVector{Float64})
    n = length(a)
    @assert length(b) == n && length(c) == n
    @assert n >= 6  # needs b[2:6] and b[n-5:n-1]

    @inbounds for i in 2:(n-1)
        num =
            (b[i+1]^2) * (a[i-1] - a[i]) -
            2.0 * b[i] * ( b[i+1]*(a[i-1]-a[i]) + b[i-1]*(a[i]-a[i+1]) ) +
            (b[i]^2)   * (a[i-1] - a[i+1]) +
            (b[i-1]^2) * (a[i] - a[i+1])

        den = (b[i-1]-b[i]) * (b[i-1]-b[i+1]) * (b[i]-b[i+1])
        c[i] = Float64(num / den)
    end

    @views begin
        AA = hchir(b[2:6], c[2:6], 5, 4)  # cubic (M=4)
        x  = Float64(b[1])
        c[1] = AA[1] + AA[2]*x + AA[3]*x^2 + AA[4]*x^3
    end

    @views begin
        AA = hchir(b[(n-5):(n-1)], c[(n-5):(n-1)], 5, 4)
        x  = Float64(b[n])
        c[n] = AA[1] + AA[2]*x + AA[3]*x^2 + AA[4]*x^3
    end

    return c
end


"""
maxpoint(b,c) -> (maxb, j)

Finds j = argmax(c) and returns maxb=b[j], j.
(Fortran maxloc(c) behavior)
"""
function maxpoint(b::AbstractVector{<:Real}, c::AbstractVector{<:Real})
    j = argmax(c)
    return (Float64(b[j]), j)
end
