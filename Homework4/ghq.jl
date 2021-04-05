using LinearAlgebra
using Plots
using DelimitedFiles

function getGaussQuadNodesWeights(alpha, beta, n)
    J = zeros((n,n))
    for i in 1:n
        if i == 1
            J[i, [1,2]] = [alpha[i], sqrt(beta[i+1])]
        elseif i == n
            J[i, [n-1, n]] = [sqrt(beta[n]), alpha[n]]
        else
            J[i, i-1:i+1] = [sqrt(beta[i]), alpha[i], sqrt(beta[i+1])]
        end
    end
    # println(J)
    F = eigen(J)
    t = F.values
    v = F.vectors
    w = beta[1] .* (v[1, :] .^ 2)
    return t, w

end

function gauss_hermite(f, n)
    alpha = zeros((n,))
    beta = convert.(Float64, collect(1:n))
    beta = (beta .- 1) ./ 2.0
    beta[1] = sqrt(pi)
    t, w = getGaussQuadNodesWeights(alpha, beta, n)
    # println(t)
    # println(w)


    # calculate the quadrature
    return sum(w .* f(t) )
end

n_end = 30
results = zeros((n_end-1, ))
mu = 0.15
sigma = 0.25
f(x) = sqrt(2) .* sigma .* 2 .* 
        (1 .+ exp.(sqrt(2) .* sigma .* x .+ mu)) .^ 0.5
exact = gauss_hermite(f, 1000)
n = 2:n_end
for (i, n) in enumerate(n)
    # get Gaussian-Hermite Quadrature Nodes and Weights
    println(n)
    results[i] = gauss_hermite(f, n)
    println(results[i])
    println()
end
@show results
@show error = abs.(results .- exact)
display(plot(collect(2:n_end), results,
            title="Gaussian-Hermite Quadrature Estimation",
            xlabel="n",
            ylabel="Estimation"))
png("GH_result.png")
display(plot(collect(n), error,
            seriestype= :scatter,
            yscale=:log10,
            ylim=[10^-15, 10^0],
            xlabel="n",
            ylabel="Error",
            title="Error of Quadrature",
            legend=false))
png("GH_error.png")

writedlm("error.csv", error)
