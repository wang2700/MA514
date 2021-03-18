using LinearAlgebra
using Plots

function spline_interpolation(f, x, xn)
    n = length(xn)
    A = zeros((n, n))
    b = zeros((n,))
    f_inter = zeros((length(x),))
    # Compute values of m
    dx = xn[2:n] .- xn[1:n-1]
    df = (f[2:n] .- f[1:n-1]) ./ dx
    A[1, 1:2] = [2, 1]
    b[1] = 3 * df[1]
    for i in range(2, stop = n - 1)
        A[i, i-1:i+1] = [dx[i], 2 * (dx[i-1] + dx[i]), 2 * dx[i-1]]
        b[i] = 3 * (dx[i] * df[i-1] + dx[i-1] * df[i])
    end
    A[n, [n - 1, n]] = [1, 2]
    b[n] = 3 * df[n-1]
    m = inv(A) * b
    # compute values of c0 to c3
    c0 = f[1:n-1]
    c1 = m[1:n-1]
    c3 = (m[2:n] .+ c1 .- 2 .* df) ./ (dx .^ 2)
    c2 = (df .- m[1:n-1]) ./ dx .- c3 .* dx
    # Calculate x
    j = 1
    for i in range(1, stop = n - 1)
        while (j <= length(x)) && (x[j] <= xn[i+1])
            f_inter[j] =
                c0[i] +
                c1[i] * (x[j] - xn[i]) +
                c2[i] * (x[j] - xn[i]) .^ 2 +
                c3[i] * (x[j] - xn[i]) .^ 3
            j = j + 1
        end
    end
    return f_inter
end

println("\n\nStart of interpolation")
println("Spline interpolation of f(x)=e^-x")
f(x) = exp.(-x)
xn = range(0, stop = 1, length = 11)
fn = f(xn)
x = range(0, stop = 1, length = 1000)
fx = f(x)
f_inter = spline_interpolation(fn, x, xn)
inter_error = abs.(f_inter - fx)
max_error = findmax(inter_error)
println(max_error)
println(x[max_error[2]])
plot(x, fx)
display(plot!(x, f_inter))

println("Spline interpolation of f(x)=x^(5/2)")
f(x) = x .^ (5 / 2)
xn = range(0, stop = 1, length = 11)
fn = f(xn)
x = range(0, stop = 1, length = 1000)
fx = f(x)
f_inter = spline_interpolation(fn, x, xn)
inter_error = abs.(f_inter - fx)
max_error = findmax(inter_error)
println(max_error)
println(x[max_error[2]])
plot(x, fx)
display(plot!(x, f_inter))

println("Spline interpolation of f(x)=x^(5/2)")
f(x) = x .^ (5 / 2)
xn = range(1, stop = 11)
xn = ((xn .- 1) ./ 10) .^ 2
fn = f(xn)
x = range(0, stop = 1, length = 1000)
fx = f(x)
f_inter = spline_interpolation(fn, x, xn)
inter_error = abs.(f_inter - fx)
max_error = findmax(inter_error)
println(max_error)
println(x[max_error[2]])
plot(x, fx)
display(plot!(x, f_inter))

println("Spline interpolation of f(x)=1/(1+x^2)")
f(x) = 1 ./ (1 .+ x .^ 2)
xn = range(-5, stop = 5, length = 25)
fn = f(xn)
x = range(0, stop = 1, length = 10000)
fx = f(x)
f_inter = spline_interpolation(fn, x, xn)
inter_error = abs.(f_inter - fx)
max_error = findmax(inter_error)
println(max_error)
println(x[max_error[2]])
plot(x, fx)
display(plot!(x, f_inter))
