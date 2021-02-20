using Pkg
Pkg.add("Plots")
using Plots
using LinearAlgebra

n_vec = [11, 101, 1001, 10001]

# compute the Chebyshev points of first kind
for n in n_vec
    k = 1:n
    x1_k = cos.((2 .* k .- 1) ./ (2 .* n) .* pi)
    plot(
        k,
        x1_k,
        xlabel = "k",
        ylabel = "Chebyshev Point Value",
        title = string("Chebyshev 1st Kind, n = ", n),
        legend = false,
        fmt = png,
    )
    png(string("1st_point_", n, ".png"))
end

# compute the Chebyshev points of second kind
for n in n_vec
    k = 0:n
    x2_k = cos.(k .* (pi / n))
    plot(
        k,
        x2_k,
        xlabel = "k",
        ylabel = "Chebyshev Point Value",
        title = string("Chebyshev 2nd Kind, n = ", n),
        legend = false,
        fmt = png,
    )
    png(string("2nd_point_", n, ".png"))
end
