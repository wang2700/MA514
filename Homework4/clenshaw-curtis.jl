using FFTW

function clenshaw_curtis(f, n)
    #generate chebyshev points
    x = cos.(pi .* collect(0:n) ./ n)
    #calculate function value at each chebyshev points
    fx = f(x) ./ (2 .* n)
    #take the fourier transform of the calculated function values
    g = real(fft(vcat(fx[1:n+1], reverse(fx[2:n]))))
    #get the chebyshev coefficients from the fft result
    a = vcat(vcat(g[1], g[2:n] .+ reverse(g[n+2:2*n])), g[n+1])
    #calcualte the weight of each chebyshev coefficient
    w = 0 .* a
    w[range(1, stop=n, step=2)] = 2 ./ 
                    (1 .- collect(range(0, stop=n, step=2)) .^ 2)
    return sum(w .* a)
end

f(x) = x .^ 2
@show clenshaw_curtis(f, 5)