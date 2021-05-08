using LinearAlgebra
using Plots

function forward_euler(f::Function, t_start::Float64, 
                    t_end::Float64, y_init, steps::Int, h::Float64)
    y = zeros(steps)
    dt = (t_end - t_start) / steps
    t = t_start
    y[1] = y_init
    for i in 1:steps - 1
        y[i + 1] = y[i] + dt .* f(t, y[i], h)
        t += dt
    end
    return y
end

function RK4(f::Function, t_start::Float64, t_end::Float64, y_init, steps::Int, h::Float64)
    y = zeros(steps)
    dt = (t_end - t_start) / steps
    t = t_start
    y[1] = y_init
    for i in 1:steps - 1
        k1 = f(t, y[i], h)
        k2 = f(t + 0.5 * dt, y[i] .+ 0.5 .* dt .* k1, h)
        k3 = f(t + 0.5 * dt, y[i] .+ 0.5 .* dt .* k2, h)
        k4 = f(t + dt, y[i] .+ dt .* k3, h)
        y[i + 1] = y[i] .+ dt ./ 6. .* (k1 + 2 .* k2 + 2 .* k3 + k4) 
        t += dt
    end
    return y
end

function f(t, y, h)
    return (sin(t + h) - 2 * sin(t) + sin(t - h)) / h^2
end

t_start = 0.
t_end = convert(Float64, pi)
println("Forward Euler")
for N in [10, 100, 1000, 10000]
    ts = range(t_start, stop=t_end, length=N)
    exact = cos.(ts)
    p = plot()
    for (i, h) in enumerate([pi / N, 0.001, 2^(-52. / 3.)])
        g_init = (sin(t_start + h) - sin(t_start - h)) / (2 * h)
        g = forward_euler(f, t_start, t_end, g_init, N, h)
        @show findmax(abs.(g - exact))
        if i == 1
            plot!(ts, g, label="N=$(N),h=pi/N")
        elseif i == 2
            plot!(ts, g, label="N=$(N),h=0.001")
        else
            plot!(ts, g, label="N=$(N),h=2^(-52/3)",
                xlabel="t",
                ylabel="g(t)")
        end
    end
    display(p)
    png("FE_N$(N).png")
end

println("RK4")
for N in [10, 100, 1000, 10000]
    ts = range(t_start, stop=t_end, length=N)
    exact = cos.(ts)
    p = plot()
    for (i, h) in enumerate([pi / N, 0.001, 2^(-52. / 3.)])
        g_init = (sin(t_start + h) - sin(t_start - h)) / (2 * h)
        g = RK4(f, t_start, t_end, g_init, N, h)
        @show findmax(abs.(g - exact))
        if i == 1
            plot!(ts, g, label="N=$(N),h=pi/N")
        elseif i == 2
            plot!(ts, g, label="N=$(N),h=0.001")
        else
            plot!(ts, g, label="N=$(N),h=2^(-52/3)",
                xlabel="t",
                ylabel="g(t)")
        end
    end
    display(p)
    png("RK_N$(N).png")
end