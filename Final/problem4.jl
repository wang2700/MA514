using LinearAlgebra
using DelimitedFiles
using Plots

"""
Compute the divided difference tables to get the Newton coefficients.
This function stores the whole table for pedagogical reasons, not
a good general purpose implementation.

The final set of coefficients to use are those on the diagonal
"""
function compute_newton_coeffs(xi, fi)
    npts = length(xi)
    T = zeros(npts, npts)
    T[:,1] = fi # the first column becomes fi.
    for j = 2:npts
        for i = j:npts
            # the value in the table is
            # made by looking left (i,j-1) and left-up (i-1,j-1)
            # divided by the region x_i - ...
            T[i,j] = (T[i,j - 1] - T[i - 1,j - 1]) / (xi[i] - xi[i - j + 1])
        end
    end
    return T
end

function compute_difference(xi, fi, pos)
    if length(pos) == 2
        return (fi[pos[2]] - fi[pos[1]]) / (xi[pos[2]] - xi[pos[1]])
    else
        n = length(pos) - 1
        f1 = compute_difference(xi, fi, pos[1:n])
        f2 = compute_difference(xi, fi, pos[2:n + 1])
        return (f2 - f1) / (xi[last(pos)] - xi[pos[1]])
    end
end

function compute_newton_coeffs_efficient(xi, fi)
    npts = length(xi)
    T = zeros(npts, )
    T[1] = fi[1]
    for i in 2:npts
        T[i] = compute_difference(xi, fi, collect(1:i))
    end
    return T
end

# f(x) = exp.(- x.^2)
# xi = range(0.0, stop=2.0, length=10)
# fi = f(xi)
# @show compute_newton_coeffs(xi, fi)
# @show compute_newton_coeffs_efficient(xi, fi)

function newton_interp(xx, T, xi)
    xv = ones(size(xx))
    f = T[1] .* xv
    for i = 2:length(T)
        xv = xv .* (xx .- xi[i - 1])
        f += T[i] * xv
    end
    return f
end

f(t) = 1 ./ (1 .+ t.^2)
r = [-5., 5.]
t = range(r[1], stop=r[2], length=500)
p = plot(t, f(t), label="Exact",
        xlabel="t",
        ylabel="f(t)")
for n in [2, 4, 6, 8]
    if (n == 2)
        ti = r
    else
        ti = range(r[1], stop=r[2], length=n)
    end
    fi = f(ti)
    coeff = compute_newton_coeffs_efficient(ti, fi)
    result = newton_interp(t, coeff, ti)
    plot!(t, result, label=string("n=", n))
end
display(p)
png("p4.png")

ti = range(r[1], stop=r[2], length=100)
fi = f(ti)
coeff = compute_newton_coeffs_efficient(ti, fi)
display(plot(ti, coeff, legend=false,
            xlabel="t",
            ylabel="Newton Coefficient"))
