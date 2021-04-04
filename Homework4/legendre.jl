using LinearAlgebra

function gauss_legendre_nodes(n)
    beta = 0.5 .* sqrt.(1 .- (2 .* collect(1:n)) .^ (-2))
    T = diagm(1 => beta, -1 => beta)
    F = eigen(T)
    x = F.values
    w = 2 .* F.vectors[1,:] .^ 2
    return x, w
end

x, w = gauss_legendre_nodes(11)
println(x)
