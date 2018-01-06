Libdl.dlopen("k_partition.so")

next_asc_k_partition(a, n, k) = ccall(:next_asc_k_partition, Int32,
    (Ptr{Int32}, Int32, Int32), pointer(a), n, k)

a = Int32[1, 1, 1, 1, 1, 1, 1, 1, 40]
n = sum(a)
k = length(a)


i = 0
@time while true
    i += 1
    println(i, ": ", a)
    if next_asc_k_partition(a, n, k) == 0
        break
    end
end


next_desc_k_partition(a, n, k) = ccall(:next_desc_k_partition, Int32,
    (Ptr{Int32}, Int32, Int32), pointer(a), n, k)

a = Int32[40, 1, 1, 1, 1, 1, 1, 1, 1]
n = sum(a)
k = length(a)

i = 0
@time while true
    i += 1
    println(i, ": ", a)
    if next_desc_k_partition(a, n, k) == 0
        break
    end
end
