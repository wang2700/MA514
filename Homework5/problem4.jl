using LinearAlgebra
using Plots

"""
This function will use Hune's method to simulate the Raptors problem
"""
function simulate_raptors_Hune(angle)
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
        dh_1, dr0_1, dr1_1, dr2_1 = compute_derivatives(angle,h,r0,r1,r2)
        dh_2, dr0_2, dr1_2, dr2_2 = compute_derivatives(angle, 
                                                        h + dt * dh_1,
                                                        r0 + dt * dr0_1,
                                                        r1 + dt * dr1_1,
                                                        r2 + dt * dr2_1)

        h += 0.5 * dt * (dh_1 + dh_2)
        r0 += 0.5 * dt * (dr0_1 + dr0_2)
        r1 += 0.5 * dt * (dr1_1 + dr1_2)
        r2 += 0.5 * dt * (dr2_1 + dr2_2)
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

N=1000
t_range = pi
angles = LinRange(-t_range, t_range, N)
tend = zeros((N,))
for (i, angle) in enumerate(angles)
    tend[i] = simulate_raptors_Hune(angle)
end
display(plot(angles, tend,
            xlabel="Angle",
            ylabel="t",
            legend=false))
png("problem4.png")