using Plots
using LinearAlgebra

# part 2
function solve_u(e::BigFloat, t::BigFloat, eps::BigFloat)
    u_n = t
    u_n_next = u_n - (u_n - e * sin(u_n) - t) * (1 - e * cos(u_n))
    while (abs(u_n - u_n_next) > eps)
        u_n = u_n_next
        u_n_next = u_n - (u_n - e * sin(u_n) - t) * (1 - e * cos(u_n))
    end
    return float(u_n_next)
end

function solve_solution(e::Float64, tspan, N)
    t = collect(LinRange(tspan[1], tspan[2], N))
    x = zeros((N,))
    y = zeros((N,))
    r = zeros((N,))
    for i in 1:N
        u = solve_u(BigFloat(e), BigFloat(t[i]), BigFloat(eps(Float64)))
        x[i] = cos(u) - e
        y[i] = sqrt(1 - e ^ 2.) * sin(u)
        r[i] = x[i] ^ 2. + y[i] ^ 2.
    end
    return x, y, r
end

tspan = (0.0, 1.0)
N = 1000
t = collect(LinRange(tspan[1], tspan[2], N))
x = zeros((N, 3))
y = zeros((N, 3))
r = zeros((N, 3))
for (i, e) in enumerate([0.3, 0.5, 0.7])
    x[:,i], y[:,i], r[:,i] = solve_solution(e, tspan, N)
end
display(plot(t, x, 
            label=["0.3" "0.5" "0.7"],
            xlabel="t",
            ylabel="x",
            title="Exact Solution of x",
            legend=:topleft))
png("exact_x.png")

display(plot(t, y, 
            label=["0.3" "0.5" "0.7"],
            xlabel="t",
            ylabel="y",
            title="Exact Solution of y",
            legend=:topleft))
png("exact_y.png")

display(plot(t, r, 
            label=["0.3" "0.5" "0.7"],
            xlabel="t",
            ylabel="r",
            title="Exact Solution of r",
            legend=:topleft))
png("exact_r.png")

# part 3

function jacobian_g(y, h)
    J = convert(Matrix{BigFloat},Matrix(1.0I, 4, 4))
    D = (y[1] ^ 2 + y[3] ^ 2) ^ 4
    J[1, 2] = -h
    J[3, 4] = -h
    J[2, 1] = -(y[3] ^ 2 - 5 * y[1] ^ 2) / D * h
    J[2, 3] = 6 * y[3] * y[1] / D * h
    J[4, 1] = 6 * y[3] * y[1] / D * h
    J[4, 3] = -(y[1] ^ 2 - 5 * y[3] ^ 2) / D * h
    return J
end


function f(t, y)
    y_out = zeros((4,))
    r = y[1] ^ 2 + y[3] ^ 2
    y_out[1] = y[2]
    y_out[2] = y[1] / r ^ 3
    y_out[3] = y[4]
    y_out[4] = y[3] / r ^ 3
    return y_out
end

function g(y_plus, y, h)
    output = convert(Vector{BigFloat}, zeros((4,)))
    r = (y_plus[1] ^ 2 + y_plus[3] ^ 2) ^ 3
    output[1] = y_plus[1] - h * y_plus[2] - y[1]
    output[2] = y_plus[2] - h * y_plus[1] / r - y[2]
    output[3] = y_plus[3] - h * y_plus[4] - y[3]
    output[4] = y_plus[4] - h * y_plus[3] / r - y[4]
    return output
end

function forward_euler(f::Function, t_start::Float64, t_end::Float64, y_init, steps::Int)
    y = zeros((length(y_init), steps))
    h = (t_end - t_start) / steps
    t = t_start
    y[:, 1] = y_init
    for i in 1:steps-1
        y[:, i + 1] = y[:, i] + h .* f(t, y[:, i])
        t += h
    end
    return y
end

function backward_euler(f::Function, t_start::Float64, t_end::Float64, y_init, steps::Int)
    y = convert(Matrix{BigFloat}, zeros((length(y_init), steps)))
    h = BigFloat((t_end - t_start) / steps)
    t = t_start
    y[:, 1] = convert(Vector{BigFloat}, y_init)
    threshold = BigFloat(eps(Float64))
    for i in 1:steps-1
        y_cur = y[:, i]
        y1 = inv(jacobian_g(y_cur, h))
        y2 = g(y_cur, y[:, i], h)
        y_diff = y1 * y2
        y_next = y_cur - y_diff
        error = abs.(y_cur .- y_next)
        while (any(x->x>threshold, error))
            y_cur = y_next
            y1 = inv(jacobian_g(y_cur, h))
            y2 = g(y_cur, y[:, i], h)
            y_diff = y1 * y2
            y_next = y_cur - y_diff
            error = abs.(y_cur .- y_next)
        end
        y[:, i+1] = y_next
    end
    return float(y)
end

e = 0.3
y_init = [1-e; 0.0; 0.0; sqrt((1 + e) / (1 - e))]
y = forward_euler(f, 0.0, 1.0, y_init, 50)
# y = backward_euler(f, 0.0, 1.0, y_init, 1000)
@show y




    