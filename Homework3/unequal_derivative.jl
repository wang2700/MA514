using Plots

function derivative_1st(f, x_vec, h)
    n = length(x)
    df = zeros((n, ))
    coeff = [-31, 63, -56, 24]
    for i in 1:n
        x = x_vec[i]
        xn = [x-2*h, x, x+h, x+1.5*h]
        df[i] = sum(f(xn) .* coeff) / (42 * h)
    end
    return df
end

function derivative_2nd(f, x_vec, h)
    n = length(x)
    df = zeros((n, ))
    # coeff = [17, -91, 154, -80]
    coeff = [24, -112, 168, -80]
    for i in 1:n
        x = x_vec[i]
        xn = [x-2*h, x, x+h, x+1.5*h]
        df[i] = sum(f(xn) .* coeff) / (42 * h^2)
    end
    return df
end

h = 1e-5
f(x) = exp.(x)
fp(x) = exp.(x)
x = range(0, stop=5, length=1000)
df = derivative_1st(f, x, h)
df2 = derivative_2nd(f, x, h)
plot(x, df, label="1st derivative approx")
plot!(x, df2, label="2nd derivative approx")
plot!(x, fp(x), label="analytical")
png("unequal_der.png")
