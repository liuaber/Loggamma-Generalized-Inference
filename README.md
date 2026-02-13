
---
title: "Generalized p-values and CIs for Gamma / Log-Gamma Quantiles"
author: ""
date: "`r format(Sys.Date(), '%Y-%m-%d')`"
output:
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
---

# Overview

This document implements simulation-based generalized p-values (and confidence limits) 
for inference in Gamma-family models using two independent sufficient statistics:

- \( S = \sum X_i \) (or \( S = \sum \log X_i \) for log-Gamma work)
- \( T = \frac{\text{geometric mean}}{\text{arithmetic mean}} \)

The generalized-inference (Gamma) component follows:

Weerahandi, S. & Gamage, J. (2016). 
*A general method of inference for two-parameter continuous distributions*. 
Communications in Statistics—Theory and Methods, 45(9), 2612–2625. 
DOI: 10.1080/03610926.2014.887109.

# Monte Carlo CDF of T

```{r}
Tcdf <- function(alp, t, n = 10, N = 500) {
  gtmp <- rgamma(n * N, shape = alp, scale = 1)
  gmat <- matrix(gtmp, N, n)
  As   <- rowMeans(gmat)
  P    <- apply(gmat, 1, prod)
  TVec <- P^(1/n) / As
  mean(TVec <= t, na.rm = TRUE)
}

GTAlpFunc <- function(alp, n, tratio, ui, Ncdf = 500) {
  Tcdf(alp, tratio, n, N = Ncdf) - ui
}
```

# Log-Gamma Quantile Generalized p-value

```{r}
pValLogGamQuan <- function(data, Q0, p, Nmc = 500, Ncdf = 500) {
  n <- length(data)
  tp <- prod(log(data))
  s  <- sum(log(data))
  amean <- s / n
  gmean <- tp^(1/n)
  t <- gmean / amean

  set.seed(round(runif(1), 5))

  U <- runif(Nmc)
  wVec <- rep(NA_real_, Nmc)

  for (i in 1:Nmc) {
    outalp <- uniroot(GTAlpFunc, c(.0001, 10),
                      n = n, tratio = t, ui = U[i], Ncdf = Ncdf)
    alphat <- outalp$root
    wVec[i] <- pgamma(s, shape = n * alphat,
                      scale = log(Q0) / qgamma(p, alphat, 1))
  }

  mean(wVec, na.rm = TRUE)
}
```

# Confidence Interval

```{r}
LogFunc <- function(x, data, p, cover, Nmc = 500, Ncdf = 500) {
  pValLogGamQuan(data, x, p, Nmc = Nmc, Ncdf = Ncdf) - cover
}

LogConLim <- function(data, p, cover = .95, Nmc = 500, Ncdf = 500) {
  temp <- uniroot(LogFunc, c(1.000001, 300),
                  data = data, p = p, cover = cover,
                  Nmc = Nmc, Ncdf = Ncdf)
  temp$root
}

CILogGamQuan <- function(data, p = .9, cover = .95,
                         Nmc = 500, Ncdf = 500) {
  LCL <- LogConLim(data, p, cover,     Nmc = Nmc, Ncdf = Ncdf)
  UCL <- LogConLim(data, p, 1 - cover, Nmc = Nmc, Ncdf = Ncdf)
  c(LCL, UCL)
}
```

# Example 1

```{r}
set.seed(123)
x <- rgamma(20, shape = 2, scale = 5)
CILogGamQuan(x, p = 0.9, cover = 0.95, Nmc = 500, Ncdf = 500)
```
