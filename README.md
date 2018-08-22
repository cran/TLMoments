TLMoments
=========

TLMoments is a set of functions which offer the functionality of calculating

-   Trimmed L-moments and

-   probability-weighted moments of given data,

-   parameters (GEV, GUM, GPD supported by now) from TL-moments,

-   quantiles of given parameters.

The calculation of TL-moments is written in C++ and therefore drastically faster than other implementations. All methods is built to be *type-conservative*, so the data's input type defines the output's type. The methods are able to process numeric vectors, matrices, lists, and dataframes (given a formula defining groups and response).

Short Introduction
==================

``` r
library(TLMoments)
```

    ## Loading required package: Rcpp

    ## Loading required package: hypergeo

    ## Loading required package: evd

    ## >> TLMoments version 0.5.995

Different data types are constructed:

``` r
xmat <- matrix(rnorm(100), nc = 4)
xvec <- xmat[, 3]
xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
xdat <- data.frame(station = rep(1:4, each = 25), hq = as.vector(xmat))
```

To calculate TL(0,1)-moments up to order 4:

``` r
TLMoments(xvec, leftrim = 0, rightrim = 1)
```

    ## $lambdas
    ##          L1          L2          L3          L4 
    ## -0.30196352  0.35745396  0.01897964  0.03025142 
    ## 
    ## $ratios
    ##          T1          T2          T3          T4 
    ##          NA -1.18376535  0.05309674  0.08463026

``` r
TLMoments(xmat, leftrim = 0, rightrim = 1)
```

    ## $lambdas
    ##           [,1]         [,2]        [,3]        [,4]
    ## L1 -0.53451579 -0.377565879 -0.30196352 -0.50307750
    ## L2  0.41147613  0.337259429  0.35745396  0.35879403
    ## L3  0.02098338  0.036463530  0.01897964 -0.05039424
    ## L4  0.05857989  0.005007545  0.03025142  0.05802465
    ## 
    ## $ratios
    ##           [,1]        [,2]        [,3]       [,4]
    ## T1          NA          NA          NA         NA
    ## T2 -0.76981099 -0.89324658 -1.18376535 -0.7131983
    ## T3  0.05099538  0.10811716  0.05309674 -0.1404545
    ## T4  0.14236523  0.01484775  0.08463026  0.1617213

``` r
TLMoments(xlist, leftrim = 0, rightrim = 1)
```

    ## $lambdas
    ## $lambdas[[1]]
    ##          L1          L2          L3          L4 
    ## -0.53451579  0.41147613  0.02098338  0.05857989 
    ## 
    ## $lambdas[[2]]
    ##           L1           L2           L3           L4 
    ## -0.377565879  0.337259429  0.036463530  0.005007545 
    ## 
    ## $lambdas[[3]]
    ##          L1          L2          L3          L4 
    ## -0.30196352  0.35745396  0.01897964  0.03025142 
    ## 
    ## $lambdas[[4]]
    ##          L1          L2          L3          L4 
    ## -0.50307750  0.35879403 -0.05039424  0.05802465 
    ## 
    ## 
    ## $ratios
    ## $ratios[[1]]
    ##          T1          T2          T3          T4 
    ##          NA -0.76981099  0.05099538  0.14236523 
    ## 
    ## $ratios[[2]]
    ##          T1          T2          T3          T4 
    ##          NA -0.89324658  0.10811716  0.01484775 
    ## 
    ## $ratios[[3]]
    ##          T1          T2          T3          T4 
    ##          NA -1.18376535  0.05309674  0.08463026 
    ## 
    ## $ratios[[4]]
    ##         T1         T2         T3         T4 
    ##         NA -0.7131983 -0.1404545  0.1617213

``` r
TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 1)
```

    ## $lambdas
    ##   station         L1        L2          L3          L4
    ## 1       1 -0.5345158 0.4114761  0.02098338 0.058579892
    ## 2       2 -0.3775659 0.3372594  0.03646353 0.005007545
    ## 3       3 -0.3019635 0.3574540  0.01897964 0.030251420
    ## 4       4 -0.5030775 0.3587940 -0.05039424 0.058024654
    ## 
    ## $ratios
    ##   station         T2          T3         T4
    ## 1       1 -0.7698110  0.05099538 0.14236523
    ## 2       2 -0.8932466  0.10811716 0.01484775
    ## 3       3 -1.1837653  0.05309674 0.08463026
    ## 4       4 -0.7131983 -0.14045451 0.16172135

As you can see, the type of the elements `lambdas` and `ratios` matches the input type.

To calculate parameters, the objects of `TLMoments` can be given to `parameters`:

``` r
tlm <- TLMoments(xvec, leftrim = 0, rightrim = 1)
parameters(tlm, "gev")
```

    ##         loc       scale       shape 
    ## -0.22837233  0.81767667  0.07383738

``` r
tlm <- TLMoments(xmat, leftrim = 0, rightrim = 1)
parameters(tlm, "gev")
```

    ##              [,1]       [,2]        [,3]       [,4]
    ## loc   -0.44813608 -0.3429497 -0.22837233 -0.2871349
    ## scale  0.94215052  0.7486603  0.81767667  0.8373964
    ## shape  0.06908072  0.1951480  0.07383738 -0.4072722

``` r
tlm <- TLMoments(xlist, leftrim = 0, rightrim = 1)
parameters(tlm, "gev")
```

    ## [[1]]
    ##        loc      scale      shape 
    ## -0.8193864  0.6902428 -0.1950720 
    ## 
    ## [[2]]
    ##         loc       scale       shape 
    ## -0.63532913  0.52899611 -0.09906923 
    ## 
    ## [[3]]
    ##        loc      scale      shape 
    ## -0.5504331  0.5982453 -0.1914612 
    ## 
    ## [[4]]
    ##        loc      scale      shape 
    ## -0.6439447  0.7007367 -0.5520610

``` r
tlm <- TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 1)
parameters(tlm, "gev")
```

    ##   station        loc     scale       shape
    ## 1       1 -0.4481361 0.9421505  0.06908072
    ## 2       2 -0.3429497 0.7486603  0.19514801
    ## 3       3 -0.2283723 0.8176767  0.07383738
    ## 4       4 -0.2871349 0.8373964 -0.40727217

To calculate quantiles, use `quantiles` after `parameters`:

``` r
tlm <- TLMoments(xvec, leftrim = 0, rightrim = 1)
quantiles(parameters(tlm, "gev"), c(.99, .999))
```

    ##     0.99    0.999 
    ## 4.250783 7.139327

``` r
tlm <- TLMoments(xmat, leftrim = 0, rightrim = 1)
quantiles(parameters(tlm, "gev"), c(.99, .999))
```

    ##           [,1]      [,2]     [,3]     [,4]
    ## 0.99  4.653664  5.235047 4.250783 1.453191
    ## 0.999 7.891586 10.588700 7.139327 1.645574

``` r
tlm <- TLMoments(xlist, leftrim = 0, rightrim = 1)
quantiles(parameters(tlm, "gev"), c(.99, .999))
```

    ## [[1]]
    ##     0.99    0.999 
    ## 1.276605 1.799341 
    ## 
    ## [[2]]
    ##     0.99    0.999 
    ## 1.319078 2.010767 
    ## 
    ## [[3]]
    ##     0.99    0.999 
    ## 1.279125 1.741558 
    ## 
    ## [[4]]
    ##      0.99     0.999 
    ## 0.5252159 0.5973432

``` r
tlm <- TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 1)
quantiles(parameters(tlm, "gev"), c(.99, .999))
```

    ##   station     0.99     0.999
    ## 1       1 4.653664  7.891586
    ## 2       2 5.235047 10.588700
    ## 3       3 4.250783  7.139327
    ## 4       4 1.453191  1.645574

The type of the outputs of `parameters` and `quantiles` are the same as the input's type as well.

Lastly, TLMoments is usable in `magrittr`-Syntax, so the cascading of functions can be written more readable as:

``` r
library(magrittr)

TLMoments(xvec, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>% 
  quantiles(c(.99, .999))
```

    ##     0.99    0.999 
    ## 4.250783 7.139327

``` r
TLMoments(xmat, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>% 
  quantiles(c(.99, .999))
```

    ##           [,1]      [,2]     [,3]     [,4]
    ## 0.99  4.653664  5.235047 4.250783 1.453191
    ## 0.999 7.891586 10.588700 7.139327 1.645574

``` r
TLMoments(xlist, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>% 
  quantiles(c(.99, .999))
```

    ## [[1]]
    ##     0.99    0.999 
    ## 1.276605 1.799341 
    ## 
    ## [[2]]
    ##     0.99    0.999 
    ## 1.319078 2.010767 
    ## 
    ## [[3]]
    ##     0.99    0.999 
    ## 1.279125 1.741558 
    ## 
    ## [[4]]
    ##      0.99     0.999 
    ## 0.5252159 0.5973432

``` r
TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 1) %>% 
  parameters("gev") %>% 
  quantiles(c(.99, .999))
```

    ##   station     0.99     0.999
    ## 1       1 4.653664  7.891586
    ## 2       2 5.235047 10.588700
    ## 3       3 4.250783  7.139327
    ## 4       4 1.453191  1.645574
