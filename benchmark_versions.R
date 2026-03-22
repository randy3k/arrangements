library(arrangements)
library(microbenchmark)

# 1. Full generation
cat("Benchmarking combinations(20, 10)\n")
print(microbenchmark(
    combinations(20, 10),
    times = 5
))

cat("\nBenchmarking permutations(10)\n")
print(microbenchmark(
    permutations(10),
    times = 5
))

# 2. Sampling (where tables help)
cat("\nBenchmarking combinations(300, 150, nsample = 100)\n")
print(microbenchmark(
    combinations(300, 150, nsample = 100),
    times = 5
))

cat("\nBenchmarking partitions(100, nsample = 100)\n")
print(microbenchmark(
    partitions(100, nsample = 100),
    times = 5
))

# 3. Nth BigZ (iterative updates)
# Halfway index for 1000 choose 500
idx <- ncombinations(1000, 500, bigz = TRUE) / 2
idx_str <- as.character(idx)
cat("\nBenchmarking combinations(1000, 500, index = <middle>) (BigZ case)\n")
print(microbenchmark(
    combinations(1000, 500, index = idx_str),
    times = 5
))

# 4. K-compositions (O(1) choose updates mentioned in commit)
cat("\nBenchmarking compositions(100, 10, nsample = 100)\n")
print(microbenchmark(
    compositions(100, 10, nsample = 100),
    times = 5
))

# 5. Multiset combinations (BigZ optimization)
cat("\nBenchmarking combinations(k = 30, freq = rep(2, 30), nsample = 100)\n")
print(microbenchmark(
    combinations(k = 30, freq = rep(2, 30), nsample = 100),
    times = 5
))
