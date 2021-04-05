using Statistics

function montecarlo(f::Function, n::Int, a::Float64, b::Float64)
    xi = (b .- a) .* rand(n) .+ a
    fi = (b .- a) .* mean(f(xi))
end


function calc_variance(f, n, iter)
    result = zeros((iter,))
    for i in 1:iter
        result[i] = montecarlo(f, n, -1., 1.)
    end
    return var(result)
end

n = 100
iter = 1000
f(x) = x .^ 2
@show montecarlo(f, n, -1., 1.)
# @show calc_variance(f, n, iter)

# f(x) = exp.(-(x .^ 2))
# @show calc_variance(f, n, iter)

# f(x) = exp.(-1. ./ (x .^ 2))
# @show calc_variance(f, n, iter)

# f(x) = exp.(x)
# @show calc_variance(f, n, iter)

# f(x) = 1 ./ (1 .+ 16 .* x .^ 2)
# @show calc_variance(f, n, iter)

# f(x) = abs.(x) .^ 3
# @show calc_variance(f, n, iter)
