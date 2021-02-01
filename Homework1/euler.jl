n = 1:20000
d = 1:5

sum = 0.0
gamma = 0.57721566490153286
for j = d
    println("\nn         gamma_n              c")
    global sum = 0.0
    println(j)
    for i = n
        global sum = sum + 1.0/i
        if i % 2000 == 0
            output = sum - log(i)
            c = (gamma - output) * (i ^ j)
            println(i, "    ", output, "   ", c,)
        end
    end
end