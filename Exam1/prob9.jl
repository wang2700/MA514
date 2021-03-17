"""Simulate a Polya Urn Process."""
function polya_var(N::Int)
  X = 1.0
  Z_total = 0.0
  Z_square_sum = 0.0
  for i=1:N
    Z = X/(i+1) # the frequency of balls
    # Z[i] is defined here.

    Z_total += Z
    Z_square_sum += Z^2

    # compute the update to X to move to the next step.
    X += rand() < Z
  end
  var = Z_square_sum / N - (Z_total / N)^2
  return var # you should also plan to change this.
end

using Random
Random.seed!(42)
var = polya_var(1000000000)
println(var)
