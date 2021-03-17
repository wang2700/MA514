using DataStructures
function sum_array(x)
    h = BinaryMinHeap(x)
    sum = 0.0
    while length(h) != 1
        min1 = pop!(h)
        min2 = pop!(h)
        push!(h, min1 + min2)
    end
    return pop!(h)
end

x = [1.5,2.3,5,3.15,21.02,2.03,45.6]
print(sum_array(x))
