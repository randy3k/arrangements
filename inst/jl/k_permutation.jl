Libdl.dlopen("k_permutation.so")

next_k_permutation(a, cycle, n, k) = ccall(:next_k_permutation, Int32,
    (Ptr{Int32}, Ptr{Int32}, Int32, Int32), pointer(a), pointer(cycle), n, k)

a = Int32[i for i in 0:9]
n = length(a)
k = 3
cycle = Int32[i for i in n:-1:(n-k+1)]

@time while true
    println(a)
    if next_k_permutation(a, cycle, n, k) == 0
        break
    end
end
