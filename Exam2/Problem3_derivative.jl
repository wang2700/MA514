using Plots
using LaTeXStrings

x = range(2.9, stop=3., length=1000)
f1(x) = -201. .* sin.(201 .* x)
f2(x) = -200. .* sin.(200 .* x)
f(x) = abs.(f1(x) - f2(x))
g(x) = abs.(cos.(200. .* x) - cos.(201. .* x))
display(plot(x, f(x),
            title="Difference between derivative of f1 and f2",
            xlabel="x",
            ylabel=L"|f_1'(x) - f_2'(x)|",
            legend=false))
png("3-Derivative_diff.png")
@show(f(3.))

display(plot(x, g(x),
            title="Difference between f1 and f2",
            xlabel="x",
            ylabel=L"|f_1(x) - f_2(x)|",
            legend=false))
png("3-Function_diff.png")
@show(g(3.))