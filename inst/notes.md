Some notes for future me.

```r
stirling2 <- function(n, k) {
    partitions(n, k, layout = "l") %>% 
      map_dbl(~ npermutations(freq = .) / prod(factorial(table(.[. > 0])))) %>%
      sum
}

bell <- function(n) {
    partitions(n, layout = "l") %>% 
      map_dbl(~ npermutations(freq = .) / prod(factorial(table(.[. > 0])))) %>%
      sum
}

ordered_bell <- function(n) {
    compositions(n, layout = "l") %>% 
        map_dbl(~ npermutations(freq = .)) %>% 
        sum
}
```
