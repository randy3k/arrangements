Libdl.dlopen("combination.so")

next_combination(a, n, k) = ccall(:next_combination, Int32,
    (Ptr{Int32}, Int32, Int32), pointer(a), n, k)

# a = Int32[i for i in 0:50]
a = Int32[0, 1, 1, 2, 3, 4]
n = length(a)
k = 5

@time while true
    println(a[1:k])
    if next_combination(a, n, k) == 0
        break
    end
end
