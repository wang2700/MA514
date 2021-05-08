using SpecialFunctions
using Plots

function C1(tf, D, R)
    x10 = 2.4
    return 2 / (x10 * besselj1(x10)) * exp(-(x10 / R)^2 * D * tf)
end
function C2(tf, D, R, Z)
    return (4 / pi) * exp(-(pi / Z)^2 * D * tf)
end
function F(tf, D, R, Z)
    c1 = C1(tf, D, R)
    c2 = C2(tf, D, R, Z)
    if c1 > 1
        return c2
    elseif c2 > 1
    return c1
  else
    return c1 * c2
    end
end
"""
Compute T_f, the final temperature,
of a simple model of cake batter cooking
for tf minutes, given
initial batter temperature - Ti
fixed oven temperature - To
diffusivity coefficient - D
cake pan radius - R
cake battery heigh - Z
"""
function FinalTemp(tf, Ti, To, D, R, Z)
    return To .+ (Ti - To) * F(tf, D, R, Z)
end

function functionD(tf, Ti, To, D, R, Z, Tf)
    return FinalTemp(tf, Ti, To, D, R, Z) - Tf
end

function bisection_D(f::Function, a::Float64, b::Float64, tol::Float64, 
                    tf, Ti, To, R, Z, Tf)
    ntol = ceil(log((b - a) / tol) / log(2))
    x = 0.0
    for i in 1:ntol
        x = 0.5 * (a + b)
        if f(tf, Ti, To, x, R, Z, Tf) < 0
            a = x
        else
            b = x
        end
    end
    return x
end

function bisection_tf(f::Function, a::Float64, b::Float64, tol::Float64, 
                    D, Ti, To, R, Z, Tf)
    ntol = ceil(log((b - a) / tol) / log(2))
    x = 0.0
    for i in 1:ntol
        x = 0.5 * (a + b)
        if f(x, Ti, To, D, R, Z, Tf) < 0
            a = x
        else
            b = x
        end
    end
    return x
end

# find D
Ti = 80.
To = 350.
R = 4
Z = 1
tf = 17
Tf = 203

tol = 0.0001
@show D = bisection_D(functionD, 0.0, 1.0, tol, tf, Ti, To, R, Z, Tf)

# estimate time
Ti = 70.
To = 375.
R = 3. 
Z = 3.
D = 0.005
Tf = 200
tol = 0.1
@show tf = bisection_D(functionD, 0.0, 200., tol, D, Ti, To, R, Z, Tf)