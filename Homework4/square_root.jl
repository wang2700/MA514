using Plots

function square_root(value::BigFloat)
    f(y) = y .^ 2 - value
    a = value
    b = 0
    x = 0
    n = 0
    while abs(f(x)) > eps(Float64)
        x = a - (a - b) / (f(a) - f(b)) * f(a)
        if f(x) > 0 && f(a) > 0
            a = x
        else
            b = x
        end
        n += 1
    end
    return x
end

# length = 100
# x_vec = range(1, stop=100, length=length)
# error_vec = zeros((length,))
# for i in 1:length
#     x = BigFloat(x_vec[i])
#     y = square_root(x)
#     error_vec[i] = abs(sqrt(x) - y)
# end
# plot(x_vec, error_vec, ylims=(0, 10 .^(-16)))
x = BigFloat(2.)
y = square_root(x)
println("False Position Estimation")
println(y)
println("Built-in function")
println(sqrt(x))
