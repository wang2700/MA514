using LinearAlgebra
using LaTeXStrings
using Plots

"""
This function will use forward Euler to simulate the Raptors problem
"""
function simulate_raptors(angle)
    vhuman = 6.0
    vraptor0 = 10.0 # the slow raptor velocity in m/s
    vraptor = 15.0 # 

    raptor_distance = 20.0

    raptor_min_distance = 0.2 # a raptor within 20 cm can attack
    tmax = 10.0 # the maximum time in seconds
    nsteps = 10000

    # initial positions
    h = [0.0,0.0]
    r0 = [1.0,0.0] * raptor_distance
    r1 = [-0.5,sqrt(3.) / 2.] * raptor_distance
    r2 = [-0.5,-sqrt(3.) / 2.] * raptor_distance

    # how much time el
    dt = tmax / nsteps
    t = 0.0

    hhist = zeros(2, nsteps + 1)
    r0hist = zeros(2, nsteps + 1)
    r1hist = zeros(2, nsteps + 2)
    r2hist = zeros(2, nsteps + 2)

    hhist[:,1] = h
    r0hist[:,1] = r0
    r1hist[:,1] = r1
    r2hist[:,1] = r2

    tend = tmax

    """
    This function will compute the derivatives of the
    positions of the human and the raptors
    """
    function compute_derivatives(angle, h, r0, r1, r2)
        dh = [cos(angle),sin(angle)] * vhuman
        dr0 = (h - r0) / norm(h - r0) * vraptor0
        dr1 = (h - r1) / norm(h - r1) * vraptor
        dr2 = (h - r2) / norm(h - r2) * vraptor
        return dh, dr0, dr1, dr2
    end

    for i = 1:nsteps
        dh, dr0, dr1, dr2 = compute_derivatives(angle, h, r0, r1, r2)
        h += dh * dt
        r0 += dr0 * dt
        r1 += dr1 * dt
        r2 += dr2 * dt
        t += dt

        hhist[:,i + 1] = h
        r0hist[:,i + 1] = r0
        r1hist[:,i + 1] = r1
        r2hist[:,i + 1] = r2

        if norm(r0 - h) <= raptor_min_distance ||
            norm(r1 - h) <= raptor_min_distance ||
            norm(r2 - h) <= raptor_min_distance

            # truncate the history
            hhist = hhist[:,1:i + 1]
            r0hist = r0hist[:,1:i + 1]
            r1hist = r1hist[:,1:i + 1]
            r2hist = r2hist[:,1:i + 1]
            tend = t
            break
        end
    end
    return tend
end

function interpolate_unif(f, r, n, pts)
    # generate the chebyshev points (2nd kind)
    k = 0:n
    x_k = cos.(k .* (pi / n))
    x_k = (r[2] + r[1]) / 2 .+ (r[2] - r[1]) / 2 .* x_k
    y = zero(pts)
    w = ones(Float64, (n + 1,))
    f_k = zeros((length(x_k),))
    for (i, x) in enumerate(x_k)
        f_k[i] = f(x)
    end

    for i in k .+ 1
        w[i] = prod(x_k[i] .- x_k[1:i - 1])
        w[i] *= prod(x_k[i] .- x_k[i + 1:n + 1])
    end
    for i in 1:size(pts, 1)
        wt = w .* (pts[i] .- x_k)
        y[i] = sum(f_k ./ wt) / sum(1.0 ./ wt)
    end
    return y

end

# generate original plot
theta_vec = range(-1., stop=1., length=1000)
f_theta = zeros((length(theta_vec),))
for (i, theta) in enumerate(theta_vec)
    f_theta[i] = simulate_raptors(theta)
end
display(plot(theta_vec, f_theta,
            xlabel=L"\theta",
            ylabel=L"f(\theta)",
            label="Exact",
            legend=:bottomright))

# chebyshev approximation
r = [-1.0, 1.0]
n = 51
approx = interpolate_unif(simulate_raptors, r, n, theta_vec)
display(plot!(theta_vec, approx,
            label="Approximation"))
png("problem2.png")