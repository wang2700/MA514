using FFTW

function clenshaw_curtis(f, n)
    x = cos.(pi .* collect(0:n) ./ n)
    fx = f(x) ./ (2 .* n)
    g = real(fft(vcat(fx[1:n+1], reverse(fx[2:n]))))
    a = vcat(vcat(g[1], g[2:n] .+ reverse(g[n+2:2*n])), g[n+1])
    w = 0 .* a
    w[range(1, stop=n, step=2)] = 2 ./ (1 .- collect(range(0, stop=n, step=2)) .^ 2)
    return sum(w .* a)
end

f(x) = x .^ 2
@show clenshaw_curtis(f, 5)