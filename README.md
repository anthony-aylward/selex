# selex
Utilities for statistical inference on SELEX sequencing data (powered by [multinomial regression](https://en.wikipedia.org/wiki/Multinomial_logistic_regression)).

## Installation and loading
```
library(devtools)
install_github("anthony-aylward/selex")
library(selex)
```

## Example
```
fit <- selex_multinom(example_counts, weights = c(8, 1, 2, 4, 8), ref = "a")
p_values <- selex_pvals(fit, n = 100, ncores = 2, timeout = 1)
```
