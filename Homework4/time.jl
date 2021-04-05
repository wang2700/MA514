using LinearAlgebra
using Plots
using DelimitedFiles

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

function selection_sort(arr)
    for i in 1:length(arr)-1
        min_idx = i
        for j in range(i+1, stop=length(arr))
            if arr[min_idx] > arr[j]
                min_idx = j
            end
        end
        temp = arr[min_idx]
        arr[min_idx] = arr[i]
        arr[i] = temp
    end
end
    
selection_sort(rand(10))

f(x) = x .^ 20
n = range(100, stop=1000, step=10)
times = zeros((length(n),))
times_sort = zeros((length(n),))
for i in 1:length(n)
    times[i] = @elapsed gauss_legendre_quad(f, n[i])
    times_sort[i] = @elapsed selection_sort(collect(1:n[i]))
end
plot(n, times,
    xlabel="n",
    ylabel="time (s)")
plot!(n, times_sort)
writedlm("time.csv", times)
writedlm("time_sort.csv", times_sort)
