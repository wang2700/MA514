using Plots
using Statistics
using CSV
using DataFrames

function derivative_1(f, x, h)
    n = length(x)
    df = zeros((n, ))
    df = (f(x .+ h) .- f(x .- h)) ./ (2 .* h)
    return df
end

function derivative_2(f, r, n, pts)
    k = 0:n
    x_k = cos.(k .* (pi / n))
    x_k = (r[2] + r[1]) / 2 .+ (r[2] - r[1]) / 2 .* x_k
    y = zeros((pts,))
    x = range(r[1], stop=r[2], length=pts)
    w = ones(Float64, (n + 1, ))
    f_k = f(x_k)
    for ind_x in 1:pts
        l_j = zeros((length(f_k),))
        for j in 1:length(k)
            for i in 1:length(k)
                product = 1
                for m in 1:length(k)
                    if (m != j && m != i)
                        product *= (x[ind_x] - x_k[m])  /
                                (x_k[j] - x_k[m])
                    end
                end
                if (i != j)
                    l_j[j] += product / (x_k[j] - x_k[i])
                end
            end
        end
        y[ind_x] = sum(f_k .* l_j)
    end
    return y
end

function eval_function_derivative_1(f, fp, f_num, startx, endx, h)
    error_vec = zeros((3, length(h)))
    for i in range(1, stop=length(h))
        x = range(startx, stop=endx, length=10000)
        df_analytical = fp(x)
        df = derivative_1(f, x, h[i])
        error = abs.(df_analytical - df)
        error_vec[:, i] = [findmax(error)[1], median(error), findmin(error)[1]]
    end
    display(plot(h, error_vec[1, :],
                title=string("strategy 1 max error of f", f_num),
                xlabel="h",
                ylabel="error"))
    png(string("f",f_num,"_d1_max.png"))
    display(plot(h, error_vec[2, :],
                title=string("strategy 1 median error of f", f_num),
                xlabel="h",
                ylabel="error"))
    png(string("f",f_num,"_d1_median.png"))
    display(plot(h, error_vec[3, :],
                title=string("strategy 1 min error of f", f_num),
                xlabel="h",
                ylabel="error"))
    png(string("f",f_num,"_d1_min.png"))
end

function eval_function_derivative_2(f, fp, f_num, startx, endx, n)
    error_vec = zeros((3, length(n)))
    x = range(startx, stop=endx, length=10000)
    df_analytical = fp(x)
    df_vec = zeros((length(n)+2, 10000))
    df_vec[1,:] = x
    df_vec[2,:] = df_analytical
    println(n)
    for i in range(1, stop=length(n))
        println(string("compute for n=", n[i]))
        df = derivative_2(f, [startx, endx], n[i], 10000)
        df_vec[i+2, :] = df
        error = abs.(df_analytical - df)
        error_vec[:, i] = [findmax(error)[1], median(error), findmin(error)[1]]
    end
    display(plot(n, error_vec[1, :],
                title=string("stratgy 2 max error of f", f_num),
                xlabel="n",
                ylabel="error"))
    png(string("f",f_num,"_d2_max.png"))
    display(plot(n, error_vec[2, :],
                title=string("stratgy 2 median error of f", f_num),
                xlabel="n",
                ylabel="error"))
    png(string("f",f_num,"_d2_median.png"))
    display(plot(n, error_vec[3, :],
                title=string("stratgy 2 min error of f", f_num),
                xlabel="n",
                ylabel="error"))
    png(string("f",f_num,"_d2_min.png"))
    writedlm(string("stratgy_2_error_f", f_num , ".csv"), error_vec, ',')
    writedlm(string("stratgy_2_df_f", f_num , ".csv"), transpose(df_vec), ',')

end

h = range(1e-10, stop=0.1, length=100)
n = 5:5:250
# n = 10
error_vec = zeros((3, length(h)))
# Evaluate derivative of f(x)=e^x
f(x) = exp.(x)
fp(x) = exp.(x)
eval_function_derivative_1(f, fp, 1, 0, 5, h)
eval_function_derivative_2(f, fp, 1, 0, 5, n)

# Evaluate derivative of f(x)=1/(1+x^2)
f(x) = 1 ./ (1 .+ x .^2)
fp(x) = -(2 .* x) ./ (1 .+ x .^ 2) .^ 2
eval_function_derivative_1(f, fp, 2, -5, 5, h)
eval_function_derivative_2(f, fp, 2, -5, 5, n)

# Evaluate derivative of f(x)=e^{3x}sin(200x^2)/(1+20x^2)
f(x) = exp.(3 .* x) .* sin.(200 .* x .^ 2) ./ (1 .+ 20 .* x .^ 2)
fp1(x) = 3 .* exp.(3 .* x) .* sin.(200 .* x.^2) ./
        (1 .+ 20 .* x .^ 2)
fp2(x) = 40 .* exp.(3 .* x) .* x .* sin.(200 .* x .^2) ./
        (1 .+ 20 .* x .^ 2).^2
fp3(x) = 400 .* exp.(3 .* x) .* x .* cos.(200 .* x .^2) ./
        (1 .+ 20 .* x .^ 2)
fp(x) = fp1(x) - fp2(x) + fp3(x)
eval_function_derivative_1(f, fp, 3, 0, 1, h)
eval_function_derivative_2(f, fp, 3, 0, 1, n)
