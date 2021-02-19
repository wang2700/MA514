using LinearAlgebra

n_vec = 5:5:25
println("n     error     cond")
for n in n_vec
  A = zeros(Float64, n+1, n+1)
  b = zeros(Float64, 1, n+1)
  for i in 0:n
    for j in 0:n
      A[i+1,j+1] = ((factorial(big(n))/(factorial(big(i))
                      *factorial(big(n-i)))) 
                  *(factorial(big(n))/(factorial(big(j))
                      *factorial(big(n-j))))
                  *((factorial(big(i+j))*factorial(big(2*n-i-j)))
                      /factorial(big(2*n+1))))
    end
    b[i+1] = 1 / (n+1)
  end
  c = b * inv(A)
  condition = cond(A)
  error = norm(c-ones(1, n+1), Inf)
  println(n, "   ", error, "   ", condition)
end
