using LinearAlgebra

function cake_derivs(x, D=0.75, delta=1.0, To=350.)
    dx = zeros((4,))
    dx[1] = D * (2 * x[2] - 2 * x[1]) / delta^2
    dx[2] = D * (x[1] - 2 * x[2] + x[3]) / delta^2
    dx[3] = D * (x[2] - 2 * x[3] + x[4]) / delta^2
    dx[4] = D * (x[3] - 2 * x[4] + To) / delta^2
    return dx
end

function evolvecake(f::Function, t_start::Float64, 
                    t_end::Float64, y_init, steps::Int, 
                    D=0.75, delta=1.0, To=350.)
    y = zeros((length(y_init), steps))
    h = (t_end - t_start) / steps
    t = t_start
    y[:, 1] = y_init
    for i in 1:steps - 1
        y[:, i + 1] = y[:, i] + h .* f(y[:, i], D, delta, To)
        t += h
    end
    return y
end

y_init = [70., 70., 70., 70.]
ts = range(0.0, stop=80., length=4800)
xhist = evolvecake(cake_derivs, 0.0, 80., y_init, 1 / 60)
p = plot()
for i = 1:4
    plot!(ts, xhist[i, :],
            xlabel="Time (Minutes)",
            ylabel="Temperature (F)", 
            label="x$(i - 1)",
            legend=:bottomright)
end
display(p)
png("problem7_3.png")

steps = round(Int, 80 / 0.77)
ts = range(0., stop=80., length=steps)
xhist = evolvecake(cake_derivs, 0.0, 80., y_init, steps)
p = plot()
for i = 1:4
    plot!(ts, xhist[i, :],
            xlabel="Time (Minutes)",
            ylabel="Temperature (F)", 
            label="x$(i - 1)",
            legend=:bottomright)
end
display(p)
png("problem7_4.png")