n_vec = [1,2,3,4,5,6,7,8,9,10]
for n in n_vec
    k = 0:n
    xi = cos.(k*pi/n)
    lambda_i = zeros(Float64, (n+1,))
    for i in k.+1
        lambda_i[i] = prod(xi[i] .- xi[1:i-1])
        lambda_i[i] *= prod(xi[i] .- xi[i+1:n+1])
    end
    println(n)
    println(lambda_i)
    println()
end
