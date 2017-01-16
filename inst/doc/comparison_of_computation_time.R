## ----echo=FALSE, message=FALSE-------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library(TLMoments)

## ------------------------------------------------------------------------
n <- c(25, 50, 100, 200, 500, 1000, 10000, 50000)
sapply(n, function(nn) {
  x <- evd::rgev(nn)
  check <- lmomco::lmoms(x, 4)$lambdas
  sapply(c("direct", "pwm", "recursive"), function(comp) {
    isTRUE(all.equal(TLMoment(x, order = 1:4, computation.method = comp), check, check.attributes = FALSE))
  })
})

## ------------------------------------------------------------------------
x <- evd::rgev(50)
microbenchmark::microbenchmark(unit = "relative", 
  TLMoments(x, max.order = 4, computation.method = "direct"), 
  TLMoments(x, max.order = 4, computation.method = "pwm"), 
  TLMoments(x, max.order = 4, computation.method = "recursive"), 
  lmomco::lmoms(x, 4), 
  Lmoments::Lmoments(x, returnobject = TRUE)
)
x <- evd::rgev(1000)
microbenchmark::microbenchmark(unit = "relative", 
  TLMoments(x, max.order = 4, computation.method = "direct"), 
  TLMoments(x, max.order = 4, computation.method = "pwm"), 
  TLMoments(x, max.order = 4, computation.method = "recursive"), 
  lmomco::lmoms(x, 4), 
  Lmoments::Lmoments(x, returnobject = TRUE)
)

## ------------------------------------------------------------------------
x <- evd::rgev(50)
microbenchmark::microbenchmark(unit = "relative", 
  TLMoment(x, order = 1:4, computation.method = "direct"), 
  TLMoment(x, order = 1:4, computation.method = "pwm"), 
  TLMoment(x, order = 1:4, computation.method = "recursive"), 
  lmom::samlmu(x, 4), 
  Lmoments::Lmoments(x, returnobject = FALSE)
)
x <- evd::rgev(1000)
microbenchmark::microbenchmark(unit = "relative", 
  TLMoment(x, order = 1:4, computation.method = "direct"), 
  TLMoment(x, order = 1:4, computation.method = "pwm"), 
  TLMoment(x, order = 1:4, computation.method = "recursive"), 
  lmom::samlmu(x, 4), 
  Lmoments::Lmoments(x, returnobject = FALSE)
)

## ------------------------------------------------------------------------
n <- c(25, 50, 100, 150, 200, 500, 1000, 10000)
sapply(n, function(nn) {
  x <- evd::rgev(nn)
  check <- lmomco::TLmoms(x, 4, leftrim = 0, rightrim = 1)$lambdas
  sapply(c("direct", "pwm", "recursive", "recurrence"), function(comp) {
    isTRUE(all.equal(TLMoment(x, order = 1:4, rightrim = 1, computation.method = comp), check, check.attributes = FALSE))
  })
})
sapply(n, function(nn) {
  x <- evd::rgev(nn)
  check <- lmomco::TLmoms(x, 4, leftrim = 2, rightrim = 4)$lambdas
  sapply(c("direct", "pwm", "recursive", "recurrence"), function(comp) {
    isTRUE(all.equal(TLMoment(x, order = 1:4, leftrim = 2, rightrim = 4, computation.method = comp), check, check.attributes = FALSE))
  })
})

## ------------------------------------------------------------------------
x <- evd::rgev(50)
microbenchmark::microbenchmark(unit = "relative", 
  TLMoments(x, max.order = 4, rightrim = 1, computation.method = "direct"), 
  TLMoments(x, max.order = 4, rightrim = 1, computation.method = "pwm"), 
  #TLMoments(x, max.order = 4, rightrim = 1, computation.method = "recursive"), 
  TLMoments(x, max.order = 4, rightrim = 1, computation.method = "recurrence"), 
  lmomco::TLmoms(x, 4, leftrim = 0, rightrim = 1)
)
x <- evd::rgev(1000)
microbenchmark::microbenchmark(unit = "relative", 
  TLMoments(x, 4, rightrim = 1, computation.method = "direct"), 
  TLMoments(x, 4, rightrim = 1, computation.method = "pwm"), 
  #TLMoments(x, 4, rightrim = 1, computation.method = "recursive"), 
  TLMoments(x, 4, rightrim = 1, computation.method = "recurrence"), 
  lmomco::TLmoms(x, 4, leftrim = 0, rightrim = 1)
)

## ------------------------------------------------------------------------
x <- evd::rgev(50)
microbenchmark::microbenchmark(unit = "relative", 
  TLMoments(x, max.order = 4, leftrim = 2, rightrim = 4, computation.method = "direct"), 
  TLMoments(x, max.order = 4, leftrim = 2, rightrim = 4, computation.method = "pwm"), 
  #TLMoments(x, max.order = 4, leftrim = 2, rightrim = 4, computation.method = "recursive"), 
  TLMoments(x, max.order = 4, leftrim = 2, rightrim = 4, computation.method = "recurrence"), 
  lmomco::TLmoms(x, 4, leftrim = 2, rightrim = 4)
)
x <- evd::rgev(1000)
microbenchmark::microbenchmark(unit = "relative", 
  TLMoments(x, max.order = 4, leftrim = 2, rightrim = 4, computation.method = "direct"), 
  TLMoments(x, max.order = 4, leftrim = 2, rightrim = 4, computation.method = "pwm"), 
  #TLMoments(x, max.order = 4, leftrim = 2, rightrim = 4, computation.method = "recursive"), 
  TLMoments(x, max.order = 4, leftrim = 2, rightrim = 4, computation.method = "recurrence"), 
  lmomco::TLmoms(x, 4, leftrim = 2, rightrim = 4)
)

