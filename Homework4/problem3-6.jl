using FFTW
using LinearAlgebra
using Statistics
using Plots
using SpecialFunctions

function clenshaw_curtis(f, n)
    #generate chebyshev points
    x = cos.(pi .* collect(0:n) ./ n)
    #calculate function value at each chebyshev points
    fx = f(x) ./ (2 .* n)
    #take the fourier transform of the calculated function values
    g = real(fft(vcat(fx[1:n+1], reverse(fx[2:n]))))
    #get the chebyshev coefficients from the fft result
    a = vcat(vcat(g[1], g[2:n] .+ reverse(g[n+2:2*n])), g[n+1])
    #calcualte the weight of each chebyshev coefficient
    w = 0 .* a
    w[range(1, stop=n+1, step=2)] = 2 ./ (1 .- collect(range(0, stop=n, step=2)) .^ 2)
    return sum(w .* a)
end

function montecarlo(f::Function, n::Int, a::Float64, b::Float64)
    xi = (b .- a) .* rand(n) .+ a
    fi = (b .- a) .* mean(f(xi))
    return fi
end

function gauss_legendre_nodes(n)
    beta = 0.5 ./ sqrt.(1 .- (2 .* collect(1:n)) .^ (-2))
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

function create_plot(f, nmax, nstep, exact, fnum)
    n_vec = range(1, stop=nmax, step=nstep)
    result = zeros((length(n_vec), 3))
    for i in 1:length(n_vec)
        result[i,1] = gauss_legendre_quad(f, n_vec[i])
        result[i,2] = clenshaw_curtis(f, n_vec[i])
        result[i,3] = montecarlo(f, n_vec[i], -1., 1.)
    end
    error = abs.(result .- exact)
    # @show error[:,3]
    display(plot(collect(n_vec), error,
        seriestype= :scatter,
        label=["Gauss-Legendre" "Clenshaw-Curtis" "Monte-Carlo"],
        yscale=:log10,
        ylim=[10^-17, 10^0],
        xlabel="n",
        ylabel="Error",
        title=string("f", fnum),
        legend=false))
    png(string("f", fnum, ".png"))
end

f(x) = x .^ 20
create_plot(f, 30, 1, 2. / 21., 1)

f(x) = exp.(-(x .^ 2))
create_plot(f, 30, 1, sqrt(pi) * erf(1), 2)

f(x) = exp.(-1. ./ (x .^ 2))
create_plot(f, 30, 1, 2. * (sqrt(pi) * (erf(1) - 1) + 1. / exp(1)), 3)

f(x) = exp.(x)
create_plot(f, 30, 1, exp(1) - 1 / exp(1), 4)

f(x) = 1 ./ (1 .+ 16 .* x .^ 2)
create_plot(f, 30, 1, 0.5 * atan(4), 5)

f(x) = abs.(x) .^ 3
create_plot(f, 30, 1, 0.5, 6)