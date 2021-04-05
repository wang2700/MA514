using LinearAlgebra

function gauss_legendre_nodes(n)
    beta = 0.5 ./ sqrt.(1 .- (2 .* collect(1:n) .^ (-2)))
    T = diagm(1 => beta, -1 => beta)
    F = eigen(T)
    x = F.values
    w = 2 .* F.vectors[1,:] .^ 2
    return x, w
end

function gauss_legendre_quad(f, n)
    x, w = gauss_legendre_nodes(n)
    return sum(f(x) .* w)
end

f(x) = x .^ 2
# x, w = gauss_legendre_nodes(11)
@show gauss_legendre_quad(f, 3)
