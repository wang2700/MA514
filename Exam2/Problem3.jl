using Plots
using Statistics

function derivative_1(f, x, h)
    df = (f(x + h) - f(x)) / h
    return df
end

function eval_function_derivative_1(f, fp, x, h)
    error_vec = zeros((length(h),))
    df = zeros((length(h),))
    for i in range(1, stop=length(h))
        df_analytical = fp(x)
        df[i] = derivative_1(f, x, h[i])
        error_vec[i] = abs(df_analytical - df[i])
    end
    return error_vec, df
end

v1 = 10000.
v2 = v1 + v1 * 0.01
f1(x) = cos.(v1 .* x)
f1p(x) = -v1 .* sin.(v1 .* x)
f2(x) = cos.(v1 .* x) + 0.001
f2p(x) = -v1 .* sin.(v1 .* x)
# f2(x) = cos.(v2 .* x)
# f2p(x) = -v2 .* sin.(v2 .* x)

h = 10 .^ range(-17, stop=-1, length=1000)
error1, df1 = eval_function_derivative_1(f1, f1p, 0, h)
plot(h, error1,
    xscale=:log10,
    xlim=[10^-17, 0.1],
    label="Function 1",
    legend=:topleft)
error2, df2 = eval_function_derivative_1(f2, f2p, 0, h)
display(plot!(h, error2,
            label="Function 2"))
plot(h, df1,
    xscale=:log10,
    xlim=[10^-17, 0.1],
    label="Function 1",
    legend=:bottomright)
display(plot!(h, df2,
    label="Function 2"))