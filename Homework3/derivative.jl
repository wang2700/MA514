using Plots
using Statistics

function derivative_1(f, x, h)
    n = length(x)
    df = zeros((n, ))
    df = (f(x .+ h) .- f(x .- h)) ./ (2 .* h)
    return df
end

function eval_function_derivative(f, fp, option, f_num, startx, endx, h)
    error_vec = zeros((3, length(h)))
    for i in range(1, stop=length(h))
        x = range(startx, stop=endx, length=10000)
        df_analytical = fp(x)
        if option == 1
            df = derivative_1(f, x, h[i])
        end
        error = abs.(df_analytical - df)
        error_vec[:, i] = [findmax(error)[1], median(error), findmin(error)[1]]
    end
    display(plot(h, error_vec[1, :],
                title=string("max error of f", f_num),
                xlabel="h",
                ylabel="error"))
    png(string("f",f_num,"_d",option,"_max.png"))
    display(plot(h, error_vec[2, :],
                title=string("median error of f", f_num),
                xlabel="h",
                ylabel="error"))
    png(string("f",f_num,"_d",option,"_median.png"))
    display(plot(h, error_vec[3, :],
                title=string("min error of f", f_num),
                xlabel="h",
                ylabel="error"))
    png(string("f",f_num,"_d",option,"_min.png"))
end

h = range(1e-10, stop=0.1, length=100)
error_vec = zeros((3, length(h)))
# Evaluate derivative of f(x)=e^x
f(x) = exp.(x)
fp(x) = exp.(x)
eval_function_derivative(f, fp, 1, 1, 0, 5, h)

# Evaluate derivative of f(x)=1/(1+x^2)
f(x) = 1 ./ (1 .+ x .^2)
fp(x) = (2 .* x) ./ (1 .+ x .^ 2) .^ 2
eval_function_derivative(f, fp, 1, 2, -5, 5, h)

# Evaluate derivative of f(x)=e^{3x}sin(200x^2)/(1+20x^2)
f(x) = exp.(3 .* x) .* sin.(200 .* x .^ 2) ./ (1 .+ 20 .* x .^ 2)
fp1(x) = 3 .* exp.(3 .* x) .* sin.(200 .* x.^2) ./
        (1 .+ 20 .* x .^ 2)
fp2(x) = 40 .* exp.(3 .* x) .* x .* sin.(200 .* x .^2) ./
        (1 .+ 20 .* x .^ 2).^2
fp3(x) = 400 .* exp.(3 .* x) .* x .* cos.(200 .* x .^2) ./
        (1 .+ 20 .* x .^ 2)
fp(x) = fp1(x) - fp2(x) + fp3(x)
eval_function_derivative(f, fp, 1, 3, 0, 1, h)
