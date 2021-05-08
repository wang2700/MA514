using LinearAlgebra

function forward_euler(f::Function, t_start::Float64, 
                    t_end::Float64, y_init, steps::Int, N)
    y = zeros((length(y_init), steps))
    h = (t_end - t_start) / steps
    t = t_start
    y[:, 1] = y_init
    for i in 1:steps - 1
        y[:, i + 1] = y[:, i] + h .* f(t, y[:, i], N)
        t += h
    end
    return y
end

function taylor(f::Function, t_start::Float64, 
                    t_end::Float64, y_init, steps::Int, N)
    y = zeros((length(y_init), steps))
    h = (t_end - t_start) / steps
    t = t_start
    y[:, 1] = y_init
    for i in 1:steps - 1
        y[:, i + 1] = y[:, i] + h .* (f(t, y[:, i], N) + 
                (h / 2) * dfdt(t, y[:, i], N))
        t += h
    end
    return y
end

function midpoint_euler(f::Function, t_start::Float64, 
                    t_end::Float64, y_init, steps::Int, N)
    y = zeros((length(y_init), steps))
    h = (t_end - t_start) / steps
    t = t_start
    y[:, 1] = y_init
    for i in 1:steps - 1
        y[:, i + 1] = y[:, i] + h .* 
            f(t + h / 2, y[:, i] + h / 2 * f(t, y[:, i], N), N)
        t += h
    end
    return y
end

function backward_euler(f::Function, t_start::Float64, 
                    t_end::Float64, y_init, steps::Int, N)
    y = convert(Matrix{BigFloat}, zeros((length(y_init), steps)))
    h = BigFloat((t_end - t_start) / steps)
    t = t_start
    y[:, 1] = convert(Vector{BigFloat}, y_init)
    threshold = BigFloat(eps(Float64))
    for i in 1:steps - 1
        # use Newton's method to solve the system of
        # non-linear equations from the Backward Euler
        # method formula.
        z = y[:, i]
        z_next = z - inv(JF(t, z, N, h)) * F(t, z, y[:, i], N, h)
        error = abs.(z .- z_next)
        while (any(x -> x > threshold, error))
            z = z_next
            z_diff = inv(JF(t, z, N, h)) * F(t, z, y[:, i], N, h)
            z_next = z - z_diff
            error = abs.(z .- z_next)
        end
        y[:, i + 1] = z_next
    end
    return convert(Matrix{Float64}, y)
end

function f(t, y, N)
    beta = 0.6
    gamma = 0.4
    y_out = zeros((3,))
    y_out[1] = - beta * y[1] * y[2] / N
    y_out[2] = beta * y[1] * y[2] / N - gamma * y[2]
    y_out[3] = gamma * y[2]
    return y_out
end

function dfdt(t, y, N)
    beta = 0.6
    gamma = 0.4
    dy = f(t, y, N)
    df = zeros((3,))
    df[1] = (- y[2] * dy[1] - y[1] * dy[2]) * beta / N 
    df[2] = (y[2] * dy[1] + y[1] * dy[2]) * beta / N - 
                gamma * dy[2]
    df[3] = gamma * dy[2]
    return df
end

function F(t, z, y, N, h)
    return z .- y - h .* f(t, z, N)
end

function JF(t, y, N, h)
    beta = 0.6
    gamma = 0.4
    J = zeros(3, 3)
    J[1, 1] = 1 + h * beta * y[2] / N
    J[1, 2] = h * beta * y[1] / N
    J[2, 1] = - h * beta * y[2] / N
    J[2, 2] = 1 - h * beta * y[1] / N + gamma
    J[3, 2] = - h * gamma
    J[3, 3] = 1
    return J
end

t_end = 25.
t_start = 0.
interval = 15. / (24. * 60.)
steps = convert(Int64, t_end / interval)
# forward Euler
println("Forward Euler")
for N in [25, 216, 237]
    y_init = [N - 1; 1; 0]
    y = forward_euler(f, t_start, t_end, y_init, steps, N)
    @show y[:, 2400]
end

# Taylor Expansion
println("Taylor Expansion")
for N in [25, 216, 237]
    y_init = [N - 1; 1; 0]
    y = taylor(f, t_start, t_end, y_init, steps, N)
    @show y[:, 2400]
end

# Midpoint Method
println("Midpoint Euler")
for N in [25, 216, 237]
    y_init = [N - 1; 1; 0]
    y = midpoint_euler(f, t_start, t_end, y_init, steps, N)
    @show y[:, 2400]
end

# Backward Method
println("Backward Euler")
for N in [25, 216, 237]
    y_init = [N - 1; 1; 0]
    y = backward_euler(f, t_start, t_end, y_init, steps, N)
    @show y[:, 2400]
end