using Plots

function half_derivative(f::Function, x, h::Float64)
    df = pi .^ -0.5 .* (f(x) ./ sqrt.(x) .+ 
            (f(x .+ h) .- f(x)) ./ h .* sqrt.(x))
    return df
end

# part 3
f(x) = 2. .* x .+ 2.
fp(x) = 4. .* sqrt.(x) ./ sqrt(pi) + 2. ./ sqrt.(pi .* x)
h = 10. ^ -6
x = collect(range(0.5, 2.5, length=100))
df = half_derivative(f, x, h)
error = abs.(df .- fp(x))
plot(x, df)
# display(plot!(x, fp(x),
#     title="Error of Approximation",
#     xlabel="x",
#     ylabel="Error",
#     legend=false))
display(plot(x, error,
    title="Error of Approximation",
    xlabel="x",
    ylabel="Error",
    legend=false))
png("4-Max_Error_3.png")
@show findmax(error)

# Part 4
h_vec = vcat(10. .^ range(-15, stop=-1, step=2), [0.5, 1.0, 2.0])
max_error = zeros((11,))
for (i, h) in enumerate(h_vec)
    df = half_derivative(f, x, h)
    error = abs.(df .- fp(x))
    println(string("h=", h))
    @show max_loc = findmax(error)
    @show x[max_loc[2]]
    max_error[i] = max_loc[1]
end
@show max_error
display(plot(h_vec, max_error,
            xscale=:log10,
            xlim=[10^-17, 2.],
            yscale=:log10,
            ylim=[10^-16, 1],
            title="Max Error Part 4"))
png("4-Max_Error_4.png")

# Part 5
f(x) = x .^ 2
fp(x) = 8. ./ (3. .* sqrt(pi)) .* x .^ 1.5
h_vec = vcat(10. .^ range(-15, stop=-1, step=2), [0.5, 1.0, 2.0])
max_error = zeros((11,))
for (i, h) in enumerate(h_vec)
    df = half_derivative(f, x, h)
    error = abs.(df .- fp(x))
    println(string("h=", h))
    @show max_loc = findmax(error)
    @show x[max_loc[2]]
    max_error[i] = max_loc[1]
end
@show max_error
display(plot(h_vec, max_error,
            xscale=:log10,
            xlim=[10^-17, 2.],
            title="Max Error Part 5"))
png("4-Max_Error_5.png")