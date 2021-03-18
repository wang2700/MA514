using Plots
using Statistics

function derivative_1(f, x, h)
    n = length(x)
    df = zeros((n-2, ))
    fx = f(x)
    df = (fx[3:n] .- fx[1:n-2]) ./ (2 .* h)
    return df, x[2:n-1]
end

function eval_function_derivative(f, fp, x, option, f_num)
    df_analytical = fp(x[2:length(x)-1])
    for i in range(1, stop=h_length)
        if option == 1
            df, dx = derivative_1(f, x, h[i])
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
                title=string("middle error of f", f_num),
                xlabel="h",
                ylabel="error"))
    png(string("f",f_num,"_d",option,"_median.png"))
    display(plot(h, error_vec[3, :],
                title=string("min error of f", f_num),
                xlabel="h",
                ylabel="error"))
    png(string("f",f_num,"_d",option,"_min.png"))
end

h_length = 1000
h = range(1e-10, stop=0.1, length=h_length)
error_vec = zeros((3, h_length))
# Evaluate derivative of f(x)=e^x
f(x) = exp.(x)
fp(x) = exp.(x)
x = range(0, stop=5, step=0.1)
eval_function_derivative(f, fp, x, 1, 1)

# Evaluate derivative of f(x)=1/(1+x^2)
f(x) = 1 ./ (1 .+ x .^2)
fp(x) = (2 .* x) ./ (1 .+ x .^ 2) .^ 2
x = range(0, stop=5, step=0.1)
eval_function_derivative(f, fp, x, 1, 2)

# Evaluate derivative of f(x)=e^{3x}sin(200x^2)/(1+20x^2)
f(x) = exp.(3 .* x) .* sin.(200 .* x .^ 2) ./ (1 .+ 20 .* x .^ 2)
fp1(x) = 3 .* exp.(3 .* x) .* sin.(200 .* x.^2) ./
        (1 .+ 20 .* x .^ 2)
fp2(x) = 40 .* exp.(3 .* x) .* x .* sin.(200 .* x .^2) ./
        (1 .+ 20 .* x .^ 2).^2
fp3(x) = 400 .* exp.(3 .* x) .* x .* cos.(200 .* x .^2) ./
        (1 .+ 20 .* x .^ 2)
fp(x) = fp1(x) - fp2(x) + fp3(x)
x = range(0, stop=5, step=0.1)
eval_function_derivative(f, fp, x, 1, 3)
