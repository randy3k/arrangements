Libdl.dlopen("multiset_combination.so")

next_multiset_combination(m, a, n, k) = ccall(:next_multiset_combination, Int32,
    (Ptr{Int32}, Ptr{Int32}, Int32, Int32), pointer(m), pointer(a), n ,k)

# a = Int32[i for i in 0:50]
a = Int32[5, 10, 10, 15, 17, 17]
m = copy(a)
n = length(a)
k = 3

@time while true
    println(a[1:k])
    if next_multiset_combination(m, a, n, k) == 0
        break
    end
end
