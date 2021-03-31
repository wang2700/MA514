using LinearAlgebra

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

function multivar_quad(f, Nx, Ny)
    # generate the nodes and weights in x and y direction
    beta = convert.(Float64, collect(1:Nx))
    beta = 1 ./ (4 .- (1 ./ (beta .- 1) .^ 2))
    beta[1] = 2.
    tx, wx = getGaussQuadNodesWeights(zeros((Nx,)), beta, Nx)

    beta = convert.(Float64, collect(1:Ny))
    beta = 1 ./ (4 .- (1 ./ (beta .- 1) .^ 2))
    beta[1] = 2.
    ty, wy = getGaussQuadNodesWeights(zeros((Ny,)), beta, Ny)

    # calculate gauss-legendre quadrature of the function
    result = 0.
    for i in 1:Nx
        for j in 1:Ny
            result += f(tx[i], ty[j]) * wx[i] * wy[j]
        end
    end
    return result
end

Nx = 3
Ny = 3

# f(x,y) = x .^ 4 + y .^ 2
# println("Function 1")
# result = multivar_quad(f, Nx, Ny)
# println(result)

f(x,y) = x .^ 4 .* y .^ 4
println("Function 2")
result = multivar_quad(f, Nx, Ny)
println(result)
#
# f(x,y) = exp.(x .^ 2 + y .^ 2)
# println("Function 3")
# result = multivar_quad(f, Nx, Ny)
# println(result)
#
# f(x,y) = (1 .- x .^ 2) + 100 .* (y .- x .^ 2) .^ 2
# println("Function 4")
# result = multivar_quad(f, Nx, Ny)
# println(result)
