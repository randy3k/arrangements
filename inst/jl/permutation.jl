Libdl.dlopen("permutation.so")

next_permutation(a, n) = ccall(:next_permutation, Int32,
    (Ptr{Int32}, Int32), pointer(a), n)

# a = Int32[i for i in 0:3]
a = [0, 0, 1, 2]
n = length(a)

i = 0
@time while true
    i += 1
    println(i, ": ", a)
    if next_permutation(a, n) == 0
        break
    end
end
