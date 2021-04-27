using LinearAlgebra
using Plots

function RK4_step(f::Function, t, y, h)
    k1 = f(t, y)
    k2 = f(t + h / 2, y .+ h ./ (2 .* k1))
    k3 = f(t + h / 2, y .+ h ./ (2 .* k2))
    k4 = f(t + h, y .+ h .* k3)
    return (k1 + 2 .* k2 + 2 .* k3 + k4) ./ 6.
end

function RK4(f::Function, t_start::Float64, t_end::Float64, y_init, steps::Int)
    y = zeros((length(y_init), steps))
    h = (t_start - t_end) / steps
    t_current = t_start
    y[:, 1] = y_init
    for i in 2:steps
        y[:, i] = RK4_step(f, t_current, y[:, i - 1], h)
        t_current += h
    end
    return y
end

function f(t, y)
    mu = 0.012277471
    mu_h = 1 - mu
    result = zeros((length(y),))
    D1 = ((y[1] + mu)^2 + y[3]^2)^1.5
    D2 = ((y[1] - mu_h)^2 + y[3]^2)^1.5
    result[1] = y[2]
    result[2] = y[1] + 2 * y[4] - mu_h * (y[1] + mu) / D1 - mu * (y[1] - mu_h) / D2 
    result[3] = y[4]
    result[4] = y[3] - 2 * y[2] - mu_h * y[3] / D1 - mu * y[3] / D2
    return result
end

y_init = zeros((4,))
y_init[1] = 0.994
y_init[2] = 0.0
y_init[3] = 0.0
y_init[4] = -2.00158510637908252240537862224
steps = 10
t_start = 0.0
t_end = 17.1
y = RK4(f, t_start, t_end, y_init, steps)
print(y)