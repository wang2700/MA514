using LinearAlgebra
using Pkg
Pkg.add("Plots")
using Plots

norm_list = zeros(Float64, 4)
for n in [1,2,3,4]
    A = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            A[i,j] = 1/(i+j)
        end
    end
    norm_list[n] = norm(A) * norm(inv(A))
end
println(norm_list)
plot(1:4, norm_list,
    xlabel = "n",
    ylabel = "Condition Number",
    legend = false,
    fmt = png)
savefig("condition_num.png")