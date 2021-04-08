using LinearAlgebra
using Plots
using LaTeXStrings

"""
This function will use forward Euler to simulate the Raptors problem
"""
function simulate_raptors(angle)
    vhuman=6.0
    vraptor0=10.0 # the slow raptor velocity in m/s
    vraptor=15.0 #

    raptor_distance = 20.0

    raptor_min_distance = 0.2 # a raptor within 20 cm can attack
    tmax=10.0 # the maximum time in seconds
    nsteps=1000

    # initial positions
    h = [0.0,0.0]
    r0 = [1.0,0.0]*raptor_distance
    r1 = [-0.5,sqrt(3.)/2.]*raptor_distance
    r2 = [-0.5,-sqrt(3.)/2.]*raptor_distance

    # how much time el
    dt = tmax/nsteps
    t = 0.0

    hhist = zeros(2,nsteps+1)
    r0hist = zeros(2,nsteps+1)
    r1hist = zeros(2,nsteps+2)
    r2hist = zeros(2,nsteps+2)

    hhist[:,1] = h
    r0hist[:,1] = r0
    r1hist[:,1] = r1
    r2hist[:,1] = r2

    tend = tmax

    """
    This function will compute the derivatives of the
    positions of the human and the raptors
    """
    function compute_derivatives(angle,h,r0,r1,r2)
        dh = [cos(angle),sin(angle)]*vhuman
        dr0 = (h-r0)/norm(h-r0)*vraptor0
        dr1 = (h-r1)/norm(h-r1)*vraptor
        dr2 = (h-r2)/norm(h-r2)*vraptor
        return dh, dr0, dr1, dr2
    end

    for i=1:nsteps
        dh, dr0, dr1, dr2 = compute_derivatives(angle,h,r0,r1,r2)
        h += dh*dt
        r0 += dr0*dt
        r1 += dr1*dt
        r2 += dr2*dt
        t += dt

        hhist[:,i+1] = h
        r0hist[:,i+1] = r0
        r1hist[:,i+1] = r1
        r2hist[:,i+1] = r2

        if norm(r0-h) <= raptor_min_distance ||
            norm(r1-h) <= raptor_min_distance ||
            norm(r2-h) <= raptor_min_distance

            # truncate the history
            hhist = hhist[:,1:i+1]
            r0hist = r0hist[:,1:i+1]
            r1hist = r1hist[:,1:i+1]
            r2hist = r2hist[:,1:i+1]
            tend = t
            break
        end
    end
    return tend
end

# Part 1 Plot f(theta)
theta_vec = range(-1., stop=1., length=1000)
f_theta = zeros((length(theta_vec),))
for (i, theta) in enumerate(theta_vec)
    f_theta[i] = simulate_raptors(theta)
end
display(plot(theta_vec, f_theta,
            xlabel= L"\theta",
            ylabel= L"f(\theta)",
            legend=false))
png("Problem2-1.png")

# Part 2 Gauss Quadrature
function gauss_legendre_nodes(n)
    beta = sqrt.((4. .- collect(1.:n) .^ -2.) .^ -1.)
    T = diagm(1 => beta, -1 => beta)
    F = eigen(T)
    x = F.values
    w = 2 .* F.vectors[1,:] .^ 2
    return x, w
end

n = 51
x, w = gauss_legendre_nodes(n)
GL_result = zeros((n,))
for i in 1:n
    GL_result[i] = simulate_raptors(x[i]) * w[i]
end
@show sum(GL_result)

# Part 3 Derivatives
theta = 0.5
h = Float64(10e-3)
# Forward Derivatives
println("Forward Derivatives")
@show (simulate_raptors(theta + h) - simulate_raptors(theta)) / h
# Backward Derivatives
println("Backward Derivatives")
@show (simulate_raptors(theta) - simulate_raptors(theta - h)) / h
# Center Derivatives
println("Center Derivatives")
@show (simulate_raptors(theta + h) - simulate_raptors(theta - h)) / (2 * h)

# Part 4 Bisection
function bisection(f::Function, a::BigFloat, b::BigFloat, tol::BigFloat)
    ntol = ceil(log((b - a) / tol) / log(2))
    x = BigFloat(0)
    for i in 1:ntol
        x = 0.5 * (a + b)
        if f(x) < 0
            a = x
        else
            b = x
        end
    end
    return x
end

a = BigFloat(0)
b = BigFloat(0.5)
f(x) = simulate_raptors(x) - 1.3625
tol = BigFloat(eps(Float64))
@show x = bisection(f, a, b, tol)
@show simulate_raptors(x)