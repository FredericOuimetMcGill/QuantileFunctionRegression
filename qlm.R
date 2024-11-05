#################
# Our functions #
#################

# Empirical quantile function, equation just above (2.1)
# Arguments:
#  p: a vector of points in (0,1) at which to compute the empirical quantile 
#  sample: the observed raw data
#  counts: NULL unless each observation is a pair (value, number of times this value was observed)
# Value:
#  values of the empirical quantile function at p[j], j = 1, ..., length(p)
qemp <- function(p, sample, counts = NULL) {
    if (!is.null(counts)) {
        if (length(sample) != length(counts)) stop("sample and counts should have the same length in function qemp()")
        if (sum(abs(na.omit(sample) - sort(na.omit(sample)))) > 10 ^ (-30)) stop("The vector of sample values should be sorted when counts are provided in function qemp()")
        if (sum(abs(na.omit(counts) - round(na.omit(counts)))) > 10 ^ (-30)) stop("The vector of count values should be integers in function qemp()")        
        tmp <- na.omit(cbind(sample, counts))
        sample <- tmp[, 1]
        counts <- tmp[, 2]
        K <- length(sample)
        N <- sum(counts)
    } else {
        sample <- sort(sample, na.last = NA) # NAs are automatically removed        
        K <- N <- length(sample)
        counts <- rep(1, N)
    }
    d <- cumsum(counts) / N
    n <- length(p)
    res <- rep(NA, n)
    for (i in 1:n) {
        res[i] <- sample[which.max(p[i] - d <= 0)]
    }
    return(res)
}

# Upper incomplete gamma function
# used in function lowergamma() below
# also used in function PhijTrunc() below
# If b < 0 and a is *not* an integer, then Gamma[a, b] is complex-valued (otherwise it is real-valued)
# but uppergamma(a, b) is not (it should!). Fortunately, this never occurs in the current research
# because uppergamma() is only used in function PhijTrunc() with b = qnorm(b) ^ 2 / 2 > 0
# or in function lowergamma() called by gammabar() with a = n - 1
uppergamma <- function(a, b) {
    if (a < 0) stop("Function uppergamma() not defined for a < 0.")
 # in Mathematica: Gamma[a, b] 
 # computes \Gamma(a, b) = \int_b^{\infty} t ^ {a - 1} * \exp(-t) dt, for (a, b)\in\mathbb{R}\times(0,\infty)
    if (any(b <= 0)) { # stop("Function uppergamma() is not defined when b <= 0.")
        res <- pracma::gammainc(b, a)["uppinc"]
        # same as gamma(a) - lowergamma(a, b)
    } else {
            res <- gamma(a) * (1 - pgamma(b, a, 1, lower = TRUE))
            # same as res <- expint::gammainc(a, b)
    }
    return(res)
}

# Lower incomplete gamma function
# used in function PhijTrunc() below
# also used in function gammabar() below
# If b < 0 and a is *not* an integer, then Gamma[a, 0, b] is complex-valued (otherwise it is real-valued)
# but lowergamma(a, b) returns NaN. Fortunately, this never occurs in the current research
# because lowergamma() is only used in function PhijTrunc() with b = qnorm(b) ^ 2 / 2 > 0
# or in function gammabar() with a = n - 1
lowergamma <- function(a, b) {
    if (a < 0) stop("Function lowergamma() not defined for a < 0.")
 # in Mathematica: Gamma[a, 0, b] or Gamma[a] - Gamma[a, b]
 # computes \gamma(a, b) = \int_0^b t ^ {a - 1} * \exp(-t) dt, for a real, b > 0
    n <- length(b)
    res <- rep(NA, n)
    if (any(b <= 0)) {
        for (i in 1:n) {
            # see https://en.wikipedia.org/wiki/Incomplete_gamma_function
            # \gamma(s, z) = \frac{z ^ s}{s} M(s, s + 1, -z) where M is Kummer's confluent hypergeometric function.
            res[i] <- b[i] ^ a * hypergeo::genhypergeo(a, a + 1, -b[i]) / a
        }
    } else {
        for (i in 1:n) {
            res[i] <- uppergamma(a, b[i])
        }
        res <- gamma(a) - res
        # same as res <- as.numeric(pracma::gammainc(b, a)["lowinc"])
    }
    return(res)
}


# Function needed by qopt() below
PhijOne <- function(j) {
    res <- if (j %% 2) 0 else factorial(j) / (2 ^ (j / 2) * factorial(j / 2))
    return(res)
}

# Function needed by qopt() below
# See Proposition A.1
PhijTrunc <- function(j, b) {
    res <- if (j %% 2) -uppergamma(j / 2 + 0.5, qnorm(b) ^ 2 / 2) else gamma(j / 2 + 0.5) + sign(b - 0.5) * lowergamma(j / 2 + 0.5, qnorm(b) ^ 2 / 2)
    return(2 ^ ((j / 2) - 1) * res / sqrt(pi))
}

# Orthogonal projection of the empirical quantile function
#  onto the space of nondecreasing polynomials of the
#  standard normal quantile function
#  See around equation (2.2)
# Arguments:
#  sample: the observed raw data
#  counts: NULL unless each observation is a pair (value, number of times this value was observed)
#  d: maximum degree of the polynomials
# Value:
#  the optimal quantile function
qopt <- function (sample, counts = NULL, d = 1) {

    if (d != 1) stop("Only the case d=1 has been implemented so far in function qopt().")

    if (!is.null(counts)) {
        if (length(sample) != length(counts)) stop("sample and counts should have the same length in function qopt()")
        if (sum(abs(na.omit(sample) - sort(na.omit(sample)))) > 10 ^ (-30)) stop("The vector of sample values should be sorted when counts are provided in function qopt()")
        if (sum(abs(na.omit(counts) - round(na.omit(counts)))) > 10 ^ (-30)) stop("The vector of count values should be integers in function qopt()")
        tmp <- na.omit(cbind(sample, counts))
        sample <- tmp[, 1]
        counts <- tmp[, 2]
        K <- length(sample)
        N <- sum(counts)
    } else {
        sample <- sort(sample, na.last = NA) # NAs are automatically removed        
        K <- N <- length(sample)
        counts <- rep(1, N)
    }
    
    n.current <- n.previous <- 0
    psi0 <- psi1 <- 0
    # See Remark A.1
    for (l in 1:K) {
        n.current <- n.current + counts[l]
        psi0 <- psi0 + sample[l] * (PhijTrunc(0, n.current / N) - PhijTrunc(0, n.previous / N))
        psi1 <- psi1 + sample[l] * (PhijTrunc(1, n.current / N) - PhijTrunc(1, n.previous / N))
        n.previous <- n.current
    }

    Phi0One <- PhijOne(0)
    Phi1One <- PhijOne(1)
    Phi2One <- PhijOne(2)
    
    a0 <- (psi1 * Phi1One - psi0 * Phi2One) / (Phi1One - Phi0One * Phi2One)
    a1 <- (psi0 * Phi1One - psi1 * Phi0One) / (Phi1One ^ 2 - Phi0One * Phi2One)

    if (a1 <= 0) stop("a1 should always be positive in function qopt()")
    
    FUN <- function(p) qnorm(p, mean = a0, sd = a1)

    environment(FUN) <- new.env()
    assign("a0", a0, environment(FUN))
    assign("a1", a1, environment(FUN))
    
    return(FUN)
}


# See Propostion 3.2
dbeta0hatML <- function(u, n, beta0, sigma2, w, muqxbar) {
    res <- dnorm(u, mean = beta0, sd = (sigma2 / n) * (1 + muqxbar ^ 2 / w))
    return(res)
}
# See Propostion 3.2
dbeta1hatML <- function(u, n, beta1, sigma2, w) {
    res <- dnorm(u, mean = beta1, sd = sigma2 / (n * w))
    return(res)
}
# See Propostion 3.2
dsigma2hatML <- function(u, n, sigma2) {
    res <- n * dchisq(n * u / sigma2, df = n - 2) / sigma2
    return(res)
}
# See Propostion 3.2
dbeta2hatML <- function(u, n, beta, sigmaqxbar, beta2) {
    res <- dexp(u - beta2, scale = beta / (n * sigmaqxbar))
    return(res)
}
# See Propostion 3.2
dbetahatML <- function(u, n, beta) {
    res <- dgamma(u, n - 1, scale = beta / n)
    return(res)    
}

# See Definition 3.7
dsigma2hat <- function(u, n, sigma2) {
    res <- (n - 2) * dchisq((n - 2) * u / sigma2, df = n - 2) / sigma2
    return(res)
}
# The regularized upper incomplete gamma function (See Proof of Proposition 3.4)
# used in function dbeta2hat()
# \bar{\Gamma}(a, b) = \frac{1}{\Gamma(a)} \int_b^{\infty} t^{a-1} e^{-t}dt = \Gamma(a, b) / \Gamma(a), a>0, b\in\mathbb{R}

Gammabar <- function(a, b) {
    # in Mathematica: GammaRegularized[a, b]
    if (a < 0) stop("Function Gammabar() not defined for a < 0.")
    if (b <= 0) {
        res <- uppergamma(a, b) / gamma(a)
        } else {
            res <- 1 - pgamma(b, a, 1, lower = TRUE)
        }
    return(res)
}
# See Proof of Proposition 3.4
dbeta2hat <- function(u, n, beta, sigmaqxbar, beta2) {
    u <- (u - beta2) / beta
    tmp <- n  * sigmaqxbar
    res <- tmp * exp(-tmp * u) * ((n - 1) / n) ^ (n - 1)
    res[u < 0] <- res[u < 0] * Gammabar(n - 1, pmax(0, -n * tmp * (u[u < 0])))
    return(res / beta)    
}
# See Definition 3.7 (or Propostion 3.4)
dbetahat <- function(u, n, beta) {
    cte <- (n - 1) / n
    res <- cte * dbetahatML(cte *  u, n, beta)
    # or equivalently
    # res <- dgamma(u / beta, n - 1, scale = 1 / (n - 1)) / beta
    return(res)
}







qlm <- function (formula, data) {

    cl <- match.call()
    attr(data, "terms") <- formula
    mf <- data
    mt <- attr(mf, "terms")

    name.qx <- as.character(formula[[3]])
    name.qy <- as.character(formula[[2]])
    qx <- data[[name.qx]]
    qy <- data[[name.qy]]
    
    # qx is a list of n explanatory quantile functions
    # qy is a list of n response quantile functions
    if (!is.list(qx)) stop("qx should be a list in function qlm()")
    if (!is.list(qy)) stop("qy should be a list in function qlm()")
    n <- length(qx)
    for (i in 1:n) {
        if (!is.function(qx[[i]])) stop(paste("Element number ", i, " in qx should be a function in function qlm()"))
        if (!is.function(qy[[i]])) stop(paste("Element number ", i, " in qy should be a function in function qlm()"))
    }
    muqx <- muqy <- sigmaqx <- sigmaqy <- rep(NA, n)
    for (i in 1:n) {
        # muqx[i] <- integrate(qx[[i]], lower = 0, upper = 1)$value
        muqx[i] <- get("a0", envir = environment(qx[[i]]))
        # muqy[i] <- integrate(qy[[i]], lower = 0, upper = 1)$value
        muqy[i] <- get("a0", envir = environment(qy[[i]]))
#        sigmaqx[i] <- sqrt(integrate(function(p) qx[[i]](p) ^ 2, lower = 0, upper = 1)$value - muqx[i] ^ 2)
        sigmaqx[i] <- get("a1", envir = environment(qx[[i]]))
#        sigmaqy[i] <- sqrt(integrate(function(p) qy[[i]](p) ^ 2, lower = 0, upper = 1)$value - muqy[i] ^ 2)
        sigmaqy[i] <- get("a1", envir = environment(qy[[i]]))
    }
    
    z <- qlm.fit(muqx, muqy, sigmaqx, sigmaqy)
    class(z) <- "qlm"
    z$call <- cl
    z$terms <- mt
    
    z$n <- n
    z
}

qlm.fit <- function (muqx, u, sigmaqx, s) {

    z <- list()

    n <- length(muqx)

    # Proposition 3.2
    ubar <- mean(u)
    sbar <- mean(s)
    muqxbar <- mean(muqx)
    sigmaqxbar <- mean(sigmaqx)
    muqx2 <- mean(muqx ^ 2)
    w <- muqx2 - muqxbar ^ 2

    z$muqx2 <- muqx2
    z$muqx <- muqx
    z$muqy <- u
    z$sigmaqx <- sigmaqx
    z$sigmaqy <- s

    # See formulas at top of Proposition 3.2
    z$muqybar <- ubar
    z$sigmaqybar <- sbar
    z$muqxbar <- muqxbar
    z$sigmaqxbar <- sigmaqxbar
    z$w <- w
    
    # Proposition 3.2
    beta1hat <- (mean(u * muqx) - ubar * muqxbar) / w # same as beta1hatML
    beta0hat <- ubar - beta1hat * muqxbar # same as beta0hatML
    # Definition 3.7 and Proposition 3.2
    sigma2hat <- sum((u - beta0hat - beta1hat * muqx) ^ 2) / (n - 2)

    mn <- min(s / sigmaqx)
    # Definition 3.7 and Proposition 3.2
    beta2hat <- (n * mn - sbar / sigmaqxbar) / (n - 1)
    betahat <- n * (sbar - mn * sigmaqxbar) / (n - 1)

    z$coefficients <- c("beta0hat" = beta0hat, "beta1hat" = beta1hat, "beta2hat" = beta2hat, "sigma2hat" = sigma2hat, "betahat" = betahat)

    # See Remark B.2
    sebeta2hat <- betahat / (sqrt(n * (n - 1)) * sigmaqxbar)

    # See Proposition 3.2 and Proposition 3.4
    sebeta0hat <- sqrt(sigma2hat * muqx2 / (n * w))
    sebeta1hat <- sqrt(sigma2hat / (n * w))
    sesigma2hat <- sqrt(2) * sigma2hat / sqrt(n - 2)
    sebetahat <- betahat / sqrt(n - 1)
    se <- c("beta0hat" = sebeta0hat, "beta1hat" = sebeta1hat, "beta2hat" = sebeta2hat, "sigma2hat" = sesigma2hat, "betahat" = sebetahat)

    z$se <- se
    
    # Section 3.2.2
    # H_0: beta_0 = 0
    # H_0: beta_1 = 0
    # H_0: beta_2 = 1
    # H_0: sigma^2 = 1
    # H_0: beta = 1
    tval <- c("beta0hat" = beta0hat / as.numeric(se["beta0hat"]), "beta1hat" = beta1hat / as.numeric(se["beta1hat"]), "beta2hat" = (beta2hat - 1) / betahat + 1 / (n * sigmaqxbar), "sigma2hat" = (n - 2) * sigma2hat, "betahat" = betahat)

    z$tval <- tval

    # Section 3.2.2
    pval <- c("beta0hat" = 2 * pt(abs(as.numeric(tval["beta0hat"])), df = n - 2, lower.tail = FALSE), "beta1hat" = 2 * pt(abs(as.numeric(tval["beta1hat"])), df = n - 2, lower.tail = FALSE), "beta2hat" = if (beta2hat + betahat / (n * sigmaqxbar) < 1) 0 else agop::ppareto2(as.numeric(tval["beta2hat"]), k = n - 1, s = (1 - 1 / n) / sigmaqxbar, lower.tail = TRUE),
              "sigma2hat" = pchisq((n - 2) * as.numeric(tval["sigma2hat"]), df = n - 2, lower.tail = TRUE), "betahat" = pgamma(as.numeric(tval["betahat"]), shape = n - 1, rate = n - 1, lower.tail = TRUE))
    
    z$pval <- pval

    # Same formulas as for the predicted values in Section 3.2.2 (but with \mu_{q_{x, n+1}} replaced with \mu_{q_{x, i}})
    # see Section 3.3
    z$fitted.values <- vector("list", n)
    for (i in 1:n) {
        muqyfitted <- beta0hat + beta1hat * muqx[i]
        sigmaqyfitted <- beta2hat * sigmaqx[i] + betahat
        z$fitted.values[[i]] <- function(x) qnorm(x, mean = muqyfitted, sd = sigmaqyfitted)
        environment(z$fitted.values[[i]]) <- new.env()
        assign("muqyfitted", muqyfitted, environment(z$fitted.values[[i]]))
        assign("sigmaqyfitted", sigmaqyfitted, environment(z$fitted.values[[i]]))
    }

    # Section 3.3
    z$residuals <- vector("list", n)
    for (i in 1:n) {
        muresid <- u[i] - beta0hat - beta1hat * muqx[i]
        sdresid <- s[i] - beta2hat * sigmaqx[i]
        z$residuals[[i]] <- function(x) qnorm(x, mean = muresid, sd = sdresid)
        environment(z$residuals[[i]]) <- new.env()
        assign("muresid", muresid, environment(z$residuals[[i]]))
        assign("sdresid", sdresid, environment(z$residuals[[i]]))
    }
    
    return(z)
}

confint.qlm <- function (object, parm, level = 0.95, ...) {
    cf <- coef(object)
    n <- object$n
    sigmaqxbar <- object$sigmaqxbar
    muqx2 <- object$muqx2
    w <- object$w
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    beta0hat <- cf["beta0hat"]
    beta1hat <- cf["beta1hat"]
    sigma2hat <- cf["sigma2hat"]
    betahat <- cf["betahat"]
    beta2hat <- cf["beta2hat"]
    ci <- list()
    alpha <- 1 - level

    se <- object$se

    # See Corollary 3.6
    cibeta0 <- beta0hat + c(-1, 1) * qt(1 - alpha / 2, n - 2) * se["beta0hat"]
    cibeta1 <- beta1hat + c(-1, 1) * qt(1 - alpha / 2, n - 2) * se["beta1hat"]
    cisigma2 <- (n - 2) * sigma2hat * (1 / c(qchisq(1 - alpha / 2, n - 2), qchisq(alpha / 2, n - 2)))
    cibeta2 <- c(beta2hat - betahat * (agop::qpareto2(1 - alpha / 2, k = n - 1, s = (1 - 1 / n) / sigmaqxbar) - 1 / (n * sigmaqxbar)), beta2hat - betahat * (agop::qpareto2(alpha / 2, k = n - 1, s = (1 - 1 / n) / sigmaqxbar) - 1 / (n * sigmaqxbar)))
    cibeta <- c(betahat / qgamma(1 - alpha / 2, n - 1, rate = n - 1), betahat / qgamma(alpha / 2, n - 1, rate = n - 1))
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    pct <- stats:::format.perc(a, 3)
    ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    ci[1, ] <- cibeta0
    ci[2, ] <- cibeta1
    ci[3, ] <- cibeta2
    ci[4, ] <- cisigma2
    ci[5, ] <- cibeta

    ci
}

print.qlm <- function(x, digits = max(3L, getOption("digits") - 3L)) {
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2L, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

plot.qlm <- function () {
}

summary.qlm <- function (object) {

    z <- object
    est <- z$coefficients
    se <- z$se
    tval <- z$tval
    pval <- z$pval
    ans <- z[c("call", "terms")]
    ans$coefficients.real <- cbind(Estimate = est[1:2], `Std. Error` = se[1:2], 
        `t value` = tval[1:2], `Pr(>|t|)` = pval[1:2])
    ans$coefficients.positive <- cbind(Estimate = est[3:5], `Std. Error` = se[3:5], 
        `stat value` = tval[3:5], `P-val (H0: par >= 1)` = pval[3:5])

    class(ans) <- "summary.qlm"
    ans

}

coef.qlm <- function (object) {
    cf <- object$coefficients
    cf
}

print.summary.qlm <- function (x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    resid <- x$residuals

    cat("\nCoefficients:\n")
    coefs.real <- x$coefficients.real
    printCoefmat(coefs.real, digits = digits, signif.stars = signif.stars, signif.legend = FALSE,
            na.print = "NA", ...)
    coefs.positive <- x$coefficients.positive
    printCoefmat(coefs.positive, digits = digits, signif.stars = TRUE, signif.legend = TRUE, P.values = TRUE, has.Pvalue = TRUE,
            na.print = "NA", ...)
    cat("\n")
    invisible(x)

    }

predict.qlm <- function(object, newdata) {

    # newdata should be a list, whose ith component is a
    # normal quantile function that needs to be created as follows:
    #   newdata[[i]] <- function(x) qnorm(x, mean = a0[i], sd = a1[i])
    #   environment(newdata[[i]]) <- new.env()
    #   assign("a0", a0[i], environment(newdata[[i]]))
    #   assign("a1", a1[i], environment(newdata[[i]]))
    
    n <- length(newdata)
    beta0hat <- object$beta0hat
    beta1hat <- object$beta1hat
    betahat <- object$betahat
    beta2hat <- object$beta2hat

    # See Section 3.2.2
    predictor <- vector("list", n)
    for (i in 1:n) {
        muqypredicted <- beta0hat + beta1hat * get("a0", envir = environment(newdata[[i]]))
        sigmaqypredicted <- beta2hat * get("a1", envir = environment(newdata[[i]])) + betahat
        predictor[[i]] <- function(x) qnorm(x, mean = muqypredicted, sd = sigmaqypredicted)
        environment(predictor[[i]]) <- new.env()
        assign("muqypredicted", muqypredicted, environment(predictor[[i]]))
        assign("sigmaqypredicted", sigmaqypredicted, environment(predictor[[i]]))
    }

    return(list(fit = predictor))
}

# Proposition 3.7
fmuQhat <- function(s, n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w) {
    res <- dnorm(s, mean = beta0hat + muqxnp1 * beta1hat, sd = sqrt(sigma2hat * (1 + (muqxnp1 - muqxbar) ^ 2 / w) / n))
    return(res)    
}

# used in function fsigmaQhat() below
# see Proposition 3.7
# regularized lower incomplete gamma function
# \bar{\gamma}(a, b) = \frac{1}{\Gamma(a)} \int_0^b t^{a-1}e^{-t}dt, a>0, b\in\mathbb{R}
gammabar <- function(a, b) {
    # in Mathematica: GammaRegularized[a, 0, b]
    if (a < 0) stop("Function gammabar() not defined for a < 0.")
    res <- lowergamma(a, b) / gamma(a)
    return(res)
}


# Proposition 3.7
fsigmaQhat <- function(t, n, betahat, beta2hat, sigmaqxnp1, sigmaqxbar) {
    tmp1 <- sigmaqxbar / sigmaqxnp1
    tmp2 <- n * tmp1 / betahat
    #part2 <- tmp2 * exp(-tmp2 * (t - beta2hat * sigmaqxnp1))
    part2 <- dexp(t - beta2hat * sigmaqxnp1, rate = tmp2)
    correction_factor <- ((n / (n - 1)) * (1 - tmp1)) ^ (n - 1)
    tmp6 <- n * (1 - tmp1) / (1 - 1 / (n * tmp1))
    upper <- tmp6 * (t - beta2hat * sigmaqxnp1) / betahat
    gamma_reg_part <- gammabar(n - 1, upper)
    res <- part2 * gamma_reg_part / correction_factor
    return(res)
}

# Proposition 3.7
fQhat <- function(s, t, n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w, betahat, beta2hat, sigmaqxnp1, sigmaqxbar) {
    tmp1 <- fmuQhat(s, n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w)
    tmp2 <- fsigmaQhat(t, n, betahat, beta2hat, sigmaqxnp1, sigmaqxbar)
    res <- tmp1 * tmp2
    return(res)
}

# Proposition 3.8
fmuE <- function(s, sigmahat, factor1i) {
    res <- dnorm(s, mean = 0, sd =  sigmahat * sqrt(factor1i))
    return(res)
}
# Proposition 3.8
fsigmaE <- function(t, n, factor2i, nu1hati, nu2hati) {
    thetahat <- 1 / (1 / nu2hati - 1 / nu1hati)
    term1 <- factor2i * dgamma(t, shape = n - 1, scale = nu2hati)
    term2 <- (1 - factor2i) * dexp(t, rate = 1 / nu1hati) * (thetahat / nu2hati) ^ (n - 2) * pgamma(t, shape = n - 2, scale = thetahat)
    res <- term1 + term2
    return(res)
}
# Proposition 3.8
fEhat <- function(s, t, n, sigmahat, factor1i, factor2i, nu1hati, nu2hati) {
    tmp1 <- fmuE(s, sigmahat, factor1i)
    tmp2 <- fsigmaE(t, n, factor2i, nu1hati, nu2hati)
    res <- tmp1 * tmp2
    return(res)
}

findkappas <- function(s, n, sigmahat, factor1i, factor2i, nu1hati, nu2hati, threshi) {
    val <- threshi / dnorm(s, mean = 0, sd =  sigmahat * sqrt(factor1i))
    upper <- 5
    kappas <- 1
    while (length(kappas) < 2) {
        kappas <- rootSolve::uniroot.all(f = function(t) return(fsigmaE(t, n, factor2i, nu1hati, nu2hati) - val), interval = c(0, upper), n = 1000)
        upper <- upper * 2
    }
    return(kappas)
}


f1s <- function(s, n, sigmahat, factor1i, factor2i, nu1hati, nu2hati, threshi) {
    lens <- length(s)
    res <- rep(NA, lens)
    for (j in 1:lens) {
        kappas <- sort(findkappas(s[j], n, sigmahat, factor1i, factor2i, nu1hati, nu2hati, threshi))
        res[j] <- (pgamma(kappas[2], n - 1, scale = nu2hati) - pgamma(kappas[1], n - 1, scale = nu2hati))
    }
    res <- dnorm(s, mean = 0, sd =  sigmahat * sqrt(factor1i)) * res
    return(res)
}

f2t <- function(t, n, nu1hati, thetahat) {
    res <- exp(-t / nu1hati) * pgamma(t, n - 2, scale = thetahat)
    return(res)
}

f2s <- function(s, n, sigmahat, factor1i, factor2i, nu1hati, nu2hati, threshi, thetahat) {
    ns <- length(s)
    res <- rep(NA, ns)
    for (j in 1:ns) {
        kappas <- sort(findkappas(s[j], n, sigmahat, factor1i, factor2i, nu1hati, nu2hati, threshi))
        res[j] <- cubature:::adaptIntegrate(f = f2t, lowerLimit = kappas[1], upperLimit = kappas[2], n = n, nu1hati = nu1hati, thetahat = thetahat, tol = 1e-04, absError = 1e-04)$integral
    }
    res <- dnorm(s, mean = 0, sd =  sigmahat * sqrt(factor1i)) * res
    return(res)
}

compute.residual.pvalue <- function(muresidi, sigmaresidi, n, sigmahat, factor1i, factor2i, nu1hati, nu2hati) {
    threshi <- fEhat(muresidi, sigmaresidi, n, sigmahat, factor1i, factor2i, nu1hati, nu2hati)
    thetahat <- 1 / (1 / nu2hati - 1 / nu1hati)
    tmpi <- (thetahat / nu2hati) ^ (n - 2) / nu1hati
    max.fsigmaEi <- abs(optim(par = 1, fn = fsigmaE, n = n, factor2i = factor2i, nu1hati = nu1hati, nu2hati = nu2hati, method = "L-BFGS-B", lower = 0, control = list(fnscale = -1))$value)
    s0i <- abs(sqrt(-2 * (sigmahat * factor1i) ^ 2 * log(sqrt(2 * pi) * sigmahat * factor1i * (threshi / max.fsigmaEi))))
    term1 <- cubature:::adaptIntegrate(f = f1s, lowerLimit = -s0i, upperLimit = s0i, n = n, sigmahat = sigmahat, factor1i = factor1i, factor2i = factor2i, nu1hati = nu1hati, nu2hati = nu2hati, threshi = threshi, tol = 1e-04, absError = 1e-04)$integral
    term2 <- cubature:::adaptIntegrate(f = f2s, lowerLimit = -s0i, upperLimit = s0i, n = n, sigmahat = sigmahat, factor1i = factor1i, factor2i = factor2i, nu1hati = nu1hati, nu2hati = nu2hati, threshi = threshi, thetahat = thetahat, tol = 1e-04, absError = 1e-04)$integral
    res <- 1 - factor2i * term1 - (1 - factor2i) * tmpi * term2
    return(res)
}

dsei <- function(s, sigmaqx, beta, i, integrate = FALSE, subdivisions = 1000L) {
    n <- length(sigmaqx)
    sum.sigmaqx <- sum(sigmaqx)
    tmp <- sigmaqx[i] / sum.sigmaqx
    nui2 <- beta * tmp / (n - 1)
    nui1 <- beta + nui2
    nui1inv <- 1 / nui1
    nui2inv <- 1 / nui2
    tmp2 <- nui2inv - nui1inv
    tmp3 <- (nui2inv / tmp2) ^ (n - 2)

    if (!integrate) {

        kappain <- sum.sigmaqx / (n * beta * sigmaqx[i])

        tmp4 <- rep(0, length(s))
        s1 <- s[s <= 1]
        s2 <- s[s > 1]
        
        somme <- rep(0, length(s1))
        for (k in 1:(n-2)) {
            prod <- 1
            for (l in 1:(k - 1)) prod <- prod * (n - l - 2) / (n * (n - 1))
            prod2 <- 1
            if (k < (n - 2)) for (i in 1:(n - k - 2)) prod2 <- prod2 * (s1 / i)
            somme <- somme + kappain ^ k * prod2 * prod
        }
        prod <- 1
        for (j in 1:n) prod <- prod * (exp(-s1 * kappain * (n - 1))) * nui2inv  / j
        somme <- somme * (n - 2) * nui2inv ^ 2 * prod
        res1 <- nui1inv * exp(-s1 * nui1inv) * (tmp3 - somme)
        
        somme <- rep(0, length(s2))
        for (k in 1:(n-2)) {
            prod <- 1
            for (l in 1:(k - 1)) prod <- prod * (n - l - 2) / (n * (n - 1))
            prod2 <- 1
            if (k < (n - 2)) for (i in 1:(n - k - 2)) prod2 <- prod2 * (s2 / i)
            somme <- somme + kappain ^ k * prod2 * prod
        }
        prod <- 1
        for (j in 1:n) prod <- prod * (exp(-s2 * kappain * (n - 1))) * nui2inv  / j
        somme <- somme * (n - 2) * nui2inv ^ 2 * prod
        res2 <- nui1inv * exp(-s2 * nui1inv) * (tmp3 - somme)
        
        tmp4[s <= 1] <- res1
        tmp4[s > 1] <- res2
        res <- tmp * dgamma(s, n - 1, scale = nui2) + (1 - tmp) * tmp4

    } else {
            myf <- function(v, s) dexp(s - v, rate = 1 / nui1) * dgamma(v, n - 2, scale = nui2)
            tmp4 <- rep(NA, length(s))
            for (j in 1:length(s)) tmp4[j] <- integrate(myf, lower = 0, upper = s[j], s = s[j], subdivisions = subdivisions)$value
            res <- tmp * dgamma(s, n - 1, scale = nui2) + (1 - tmp) * tmp4
        }
    return(res)
}
