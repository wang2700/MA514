k = 1:7
sum1 = 0.0
sum2 = 0.0

println("k   error1                 error2")
for i = k
    n = 1:10^i
    global sum1 = 0.0
    global sum2 = 0.0
    for j = n
        global sum1 += 1.0/j - 1.0/(j+1)
        global sum2 += 1.0/(j * (j + 1))
    end
    actual = 1 - 1.0/(10^i + 1)
    println(i, "   ", abs(actual-sum1), "   ", abs(actual-sum2))
end
