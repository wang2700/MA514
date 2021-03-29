using LinearAlgebra
using Plots

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

n_end = 10
results = zeros((n_end-4, ))
for n = 5:n_end
    # get Gaussian-Hermite Quadrature Nodes and Weights
    println(n)
    alpha = zeros((n,))
    beta = convert.(Float64, collect(1:n))
    beta = (beta .- 1) ./ 2.0
    beta[1] = sqrt(pi)
    t, w = getGaussQuadNodesWeights(alpha, beta, n)
    # println(t)
    # println(w)


    # calculate the quadrature
    mu = 0.15
    sigma = 0.25
    f(x) = 2 .* (1 .+ exp.(sqrt(2) .* sigma .* x .+ mu)) .^ 0.5
    results[n-4] = sum(w .* f(t))
    println(results[n-4])
    println()
end
display(plot(collect(5:n_end), results,
            title="Gaussian-Hermite Quadrature Estimation",
            xlabel="n",
            ylabel="Estimation"))
png("GH_result.png")
