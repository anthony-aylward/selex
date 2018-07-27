# selex
Utilities for statistical inference on SELEX sequencing data

## Installation and loading
```
library(devtools)
install_github("anthony-aylward/selex")
library(selex)
```

## Example
```
fit <- selex_multinom(example_counts, weights = c(8, 1, 2, 4, 8), ref = "a")
p_values <- selex_pvals(fit)
```
