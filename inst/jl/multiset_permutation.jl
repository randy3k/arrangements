Libdl.dlopen("multiset_permutation.so")

next_multiset_permutation(a, n, k) = ccall(:next_multiset_permutation, Int32,
    (Ptr{Int32}, Int32, Int32), pointer(a), n, k)

a = Int32[i for i in 0:9]
n = length(a)
k = length(a)

@time while true
    # println(a)
    if next_k_permutation(a, n, k) == 0
        break
    end
end
