# using Pkg
# Pkg.add("Plots")
using Plots
using LinearAlgebra

function interpolate_unif(f, r, n, pts)
    # generate the chebyshev points (2nd kind)
    k = 0:n
    x_k = cos.(k .* (pi / n))
    x_k = (r[2] + r[1]) / 2 .+ (r[2] - r[1]) / 2 .* x_k
    y = zero(pts)
    w = ones(Float64, (n + 1,))
    f_k = f(x_k)
    for i in k .+ 1
        w[i] = prod(x_k[i] .- x_k[1:i-1])
        w[i] *= prod(x_k[i] .- x_k[i+1:n+1])
    end
    for i in 1:size(pts, 1)
        wt = w .* (pts[i] .- x_k)
        y[i] = sum(f_k ./ wt) / sum(1.0 ./ wt)
    end
    return y

end

# f(x) = exp.(-x)
# plot_title = "Chebyshev Approximation of e^x"
# file_name = "chebyshev_approx_ex.png"
f(x) = 1 ./ (x .^ 2 .+ 5)
plot_title = "Chebyshev Approximation of 1/(x^2+5)"
file_name = "chebyshev_approx_poly.png"
pts = range(1, stop = 5, length = 10000)
n_vec = [5, 8, 11, 30, 51]
plot(pts, f(pts), label="function", title=plot_title)
for n in n_vec
    y = interpolate_unif(f, [1, 5], n, pts)
    display(plot!(pts,
                y,
                label=string("n=", n)))
end
png(string("Homework2/", file_name))
