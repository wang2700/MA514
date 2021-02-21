using PyPlot
using LinearAlgebra

f(x,y) = 1.0 ./ (x.^2 .+ y.^2 .+ 5)
r = [-5,5]
n = 100
x_vec = range(r[1], stop=r[2], length=n)
y_vec = range(r[1], stop=r[2], length=n)

x_grid = repeat(x_vec',n,1)
y_grid = repeat(y_vec,1,n)

z = zeros(n,n)

z_func = f(x_grid, y_grid)

fig = figure("function_surfactplot", figsize=(10,10))
plot_surface(
    x_grid,
    y_grid,
    z_func,
    rstride = 2,
    edgecolors = "k",
    cstride = 2,
    cmap = ColorMap("gray"),
    alpha = 0.8,
    linewidth = 0.25,
)
xlabel("X")
ylabel("Y")
PyPlot.title("Surface Plot of Function")
show()

# generate Chebyshev points
println("Generate Chebyshev Points")
k = 1:5
c_pts = cos.((2 .* k .- 1) ./ (2 .* 5) .* pi)
c_pts = (r[2] + r[1]) / 2 .+ (r[2] - r[1]) / 2 .* c_pts

for ind_x in 1:n
    for ind_y in 1:n
        x = x_vec[ind_x]
        y = y_vec[ind_y]
        for i = 1:5
            for j = 1:5
                Lix = 1
                Liy = 1
                for s = 1:5
                    if (i != s)
                        Lix *= (x-c_pts[s])/(c_pts[i]-c_pts[s])
                    end
                    if (j != s)
                        Liy *= (y-c_pts[s])/(c_pts[j]-c_pts[s])
                    end
                end
                L = Lix * Liy
                z[ind_x, ind_y] += f(c_pts[i], c_pts[j]) * L
            end
        end
    end
end

fig = figure("interpolation_surfactplot", figsize=(10,10))
plot_surface(
    x_grid,
    y_grid,
    z,
    rstride = 2,
    edgecolors = "k",
    cstride = 2,
    cmap = ColorMap("gray"),
    alpha = 0.8,
    linewidth = 0.25,
)
xlabel("X")
ylabel("Y")
PyPlot.title("Surface Plot of Interpolant")
show()

fig = figure("interpolation_contour", figsize=(10,10))
cp =contour(
    x_grid,
    y_grid,
    z,
    color="black",
    linewidth=2.0,
)
xlabel("X")
ylabel("Y")
PyPlot.title("Contour Plot of Interpolant")
show()
