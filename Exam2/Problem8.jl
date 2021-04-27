using Plots
using LinearAlgebra

function gauss_legendre_nodes(n)
    beta = sqrt.((4. .- collect(1.:n) .^ -2.) .^ -1.)
    T = diagm(1 => beta, -1 => beta)
    F = eigen(T)
    x = F.values
    w = 2 .* F.vectors[1,:] .^ 2
    return x, w
end

n_vec = 5:100
max_x = zeros((length(n_vec),))
max_w = zeros((length(n_vec),))
for (i, n) in enumerate(n_vec)
    x, w = gauss_legendre_nodes(n)
    max_x[i] = findmax(x)[1]
    max_w[i] = findmax(w)[1]
end
display(plot(n_vec, max_x,
            ylabel="Max x",
            xlabel="n",
            legend=false))
png("8-Max_x.png")
display(plot(n_vec, max_w,
            ylabel="Max w",
            xlabel="n",
            legend=false))
png("8-Max_w.png")


