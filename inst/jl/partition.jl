Libdl.dlopen("partition.so")

next_asc_partition(a, k) = ccall(:next_asc_partition, Int32,
    (Ptr{Int32}, Ptr{Int32}), pointer(a), k)

n = 6
# a = vcat(Int32[0, n], Int32[0 for _ in 1:n])
a = Int32[1 for _ in 1:n]
k = Ref{Int32}(n-1)

i = 0
@time while true
    i += 1
    println(i, ": ", a[1:(k[]+1)])
    if next_asc_partition(a, k) == 0
        break
    end
end


next_desc_partition(a, h, k) = ccall(:next_desc_partition, Int32,
    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}), pointer(a), h, k)

n = 6
# a = vcat(Int32[0, n], Int32[0 for _ in 1:n])
a = Int32[1 for _ in 1:n]
a[1] = Int32(n)
h = Ref{Int32}(0)
k = Ref{Int32}(1)

i = 0
@time while true
    i += 1
    println(i, ": ", a[1:(k[])])
    if next_desc_partition(a, h, k) == 0
        break
    end
end
