using LinearAlgebra
using Plots
using DelimitedFiles

function RK4(f::Function, t_start::Float64, t_end::Float64, y_init, steps::Int)
    y = zeros((length(y_init), steps))
    h = (t_end - t_start) / steps
    t = t_start
    y[:, 1] = y_init
    for i in 1:steps-1
        k1 = f(t, y[:, i])
        k2 = f(t + 0.5 * h, y[:, i] .+ 0.5 .* h .* k1)
        k3 = f(t + 0.5 * h, y[:, i] .+ 0.5 .* h .* k2)
        k4 = f(t + h, y[:, i] .+ h .* k3)
        y[:, i + 1] = y[:, i] .+ h ./ 6. .* (k1 + 2 .* k2 + 2 .* k3 + k4) 
        t += h
    end
    return y
end

function astronomy(t, y)
    mu = 0.012277471
    mu_h = 1 - mu
    result = zeros((length(y),))
    D1 = ((y[1] + mu) ^ 2 + y[3] ^ 2) ^ 1.5
    D2 = ((y[1] - mu_h) ^ 2 + y[3] ^ 2) ^ 1.5
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
t_start = 0.0
t_end = 17.1
steps = 
for steps in [1000, 5000, 10000, 50000, 100000]
    y = RK4(astronomy, t_start, t_end, y_init, steps)
    display(plot(y[1,:], y[3,:],
                xlabel="u_1",
                ylabel="u_2",
                title=string("steps = ", steps)))
    png(string("problem2_", steps, ".png"))
end
