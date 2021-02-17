# a = 1.23
# b = 1.03
# c = 0.63

x = rand(Float64, 3) * 1e200
x = sort(x, rev=true)
a = x[1]
b = x[2]
c = x[3]

@show x
println((a+b)+c)
println(a+(b+c))
println((a+c)+b)