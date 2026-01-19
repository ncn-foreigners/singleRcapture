## ----plot-parameters, include=FALSE-------------------------------------------
par(mar = c(2.5, 8.5, 4.1, 2.5), cex.main = .7, cex.lab = .6)

## ----eval=FALSE---------------------------------------------------------------
# ?ztpoisson

## ----installation, eval=FALSE-------------------------------------------------
# install.packages("singleRcapture")

## ----loading-package, warning=FALSE-------------------------------------------
library(singleRcapture)

## ----loading_singleRcapture---------------------------------------------------
data(netherlandsimmigrant)
head(netherlandsimmigrant)

## ----netherlandsimmigrant_summary---------------------------------------------
summary(netherlandsimmigrant)

## ----table_netherlands--------------------------------------------------------
table(netherlandsimmigrant$capture)

## ----basic_model_netherlands--------------------------------------------------
basicModel <- estimatePopsize(
  formula = capture ~ gender + age + nation,
  model   = ztpoisson(),
  data    = netherlandsimmigrant,
  controlMethod = controlMethod(silent = TRUE)
)
summary(basicModel)

## ----model_inflated_netherlands, warning=FALSE, cache=T-----------------------
set.seed(123456)
modelInflated <- estimatePopsize(
    formula = capture ~ nation,
    model   = oiztgeom(omegaLink = "cloglog"),
    data    = netherlandsimmigrant,
    controlModel = controlModel(
        omegaFormula = ~ gender + age
    ),
    popVar = "bootstrap",
    controlPopVar = controlPopVar(bootType = "semiparametric",
                                  B = 50)
)
summary(modelInflated)

## ----extract_pop_size---------------------------------------------------------
popSizeEst(basicModel)    # alternative: basicModel$populationSize
popSizeEst(modelInflated) # alternative: modelInflated$populationSize

## ----message=FALSE------------------------------------------------------------
library(lmtest)

## -----------------------------------------------------------------------------
lrtest(basicModel, modelInflated, 
       name = function(x) {
    if (family(x)$family == "ztpoisson")
        "Basic model"
    else "Inflated model"
})

## ----marginal_freq_basic_model------------------------------------------------
margFreq <- marginalFreq(basicModel)
summary(margFreq, df = 1, dropl5 = "group")

## ----marginal_freq_inflated_model---------------------------------------------
margFreq_inf <- marginalFreq(modelInflated)
summary(margFreq_inf, df = 1, dropl5 = "group")

## ----rootogram, fig.dim=c(5,5), fig.pos="ht", fig.show='hold', fig.cap="Rootograms for ztpoisson (left) and oiztgeom (right) models"----
plot(   basicModel, plotType = "rootogram", main = "ZT Poisson model")
plot(modelInflated, plotType = "rootogram", main = "OI ZT Geometric model")

## ----dfbeta_basic_model-------------------------------------------------------
dfb <- dfbeta(basicModel)
round(t(apply(dfb, 2, quantile)*100), 4)

## ----df_beta_inflated_model---------------------------------------------------
dfi <- dfbeta(modelInflated)
round(t(apply(dfi, 2, quantile)*100), 4)

## ----dfpopsize_basic_model----------------------------------------------------
dfb_pop <- dfpopsize(basicModel, dfbeta = dfb)
dfi_pop <- dfpopsize(modelInflated, dfbeta = dfi)
summary(dfb_pop)
summary(dfi_pop)

## ----dfpopsize_plot, fig.dim=c(5,5), fig.pos="ht", fig.show='hold', fig.cap="Results for ztpoisson (left) and oiztgeom (right) model"----
plot(basicModel, plotType = "dfpopContr", 
     dfpop = dfb_pop, xlim = c(-4500, 150))
plot(modelInflated, plotType = "dfpopContr", 
     dfpop = dfi_pop, xlim = c(-4500, 150))

## ----eval=FALSE---------------------------------------------------------------
# ?plot.singleRStaticCountData

## ----strata-------------------------------------------------------------------
popSizestrata <- stratifyPopsize(basicModel)
cols <- c("name", "Observed", "Estimated", "logNormalLowerBound", 
          "logNormalUpperBound")
popSizestrata_report <- popSizestrata[, cols]
cols_custom <- c("Name", "Obs", "Estimated", "LowerBound", "UpperBound")
names(popSizestrata_report) <- cols_custom
popSizestrata_report

## ----strata-inflated----------------------------------------------------------
popSizestrata_inflated <- stratifyPopsize(modelInflated)
popSizestrata_inflated_report <- popSizestrata_inflated[, cols]
names(popSizestrata_inflated_report) <- cols_custom
popSizestrata_inflated_report

## ----custom-strata, message=F, warning=F--------------------------------------
library(sandwich)
popSizestrataCustom <- stratifyPopsize(
  object  = basicModel,
  strata = ~ gender + age, 
  alpha   = rep(c(0.1, 0.05), each=2), 
  cov     = vcovHC(basicModel, type = "HC4")
)

popSizestrataCustom_report <- popSizestrataCustom[, c(cols, "confLevel")]
names(popSizestrataCustom_report) <- c(cols_custom, "alpha")
popSizestrataCustom_report

## ----eval=FALSE---------------------------------------------------------------
# list(
#   "Stratum 1" = netherlandsimmigrant$gender == "male"   &
#     netherlandsimmigrant$nation == "Suriname",
#   "Stratum 2" = netherlandsimmigrant$gender == "female" &
#     netherlandsimmigrant$nation == "North Africa"
# )

## ----strata_plot, fig.dim=c(6,6), fig.show='hold', fig.pos = "ht", fig.cap="Population size by covariates for ztpoisson (left) and oiztgeom (right) model"----
plot(basicModel, plotType = "strata")
plot(modelInflated, plotType = "strata")

## ----popsize_extract----------------------------------------------------------
(popEst <- popSizeEst(basicModel))

## ----extract-coefs------------------------------------------------------------
coef(summary(basicModel))

## ----simulate-method----------------------------------------------------------
set.seed(1234567890)
N <- 10000
gender <- rbinom(N, 1, 0.2)
eta <- -1 + 0.5*gender
counts <- simulate(ztpoisson(), eta = cbind(eta), seed = 1)
summary(data.frame(gender, eta, counts))

## ----example-code-farm, eval=FALSE--------------------------------------------
# estimatePopsize(
#   TOTAL_SUB ~ .,
#   data = farmsubmission,
#   model = ztoigeom(),
#   controlModel(
#     omegaFormula = ~ 1 + log_size + C_TYPE
#   )
# )

## ----estimatePopsizeFit-1-----------------------------------------------------
X <- matrix(data = 0, nrow = 2 * NROW(farmsubmission), ncol = 7)

## ----estimatePopsizeFit-2-----------------------------------------------------
X[1:NROW(farmsubmission), 1:4] <- model.matrix(
  ~ 1 + log_size + log_distance + C_TYPE, 
  farmsubmission
)
X[-(1:NROW(farmsubmission)), 5:7] <- model.matrix(
  ~ 1 + log_distance + C_TYPE, 
  farmsubmission
)
attr(X, "hwm") <- c(4, 3)

## ----estimatePopsizeFit-3-----------------------------------------------------
start <- glm.fit(
  y = farmsubmission$TOTAL_SUB, 
  x = X[1:NROW(farmsubmission), 1:4], 
  family = poisson()
)$coefficients
start

## ----estimatePopsizeFit-4-----------------------------------------------------
res <- estimatePopsizeFit(
  y            = farmsubmission$TOTAL_SUB, 
  X            = X, 
  method       = "IRLS", 
  priorWeights = 1, 
  family       = ztoigeom(), 
  control      = controlMethod(silent = TRUE), 
  coefStart    = c(start, 0, 0, 0),
  etaStart     = matrix(X %*% c(start, 0, 0, 0), ncol = 2),
  offset       = cbind(rep(0, NROW(farmsubmission)), 
                       rep(0, NROW(farmsubmission)))
)

## ----estimatePopsizeFit-5-----------------------------------------------------
ll <- ztoigeom()$makeMinusLogLike(y = farmsubmission$TOTAL_SUB, X = X)

## ----estimatePopsizeFit-6-----------------------------------------------------
res2 <- estimatePopsizeFit(
  y = farmsubmission$TOTAL_SUB, 
  X = X, 
  method = "optim", 
  priorWeights = 1, 
  family = ztoigeom(), 
  coefStart = c(start, 0, 0, 0),
  control = controlMethod(silent = TRUE, maxiter = 10000),
  offset = cbind(rep(0, NROW(farmsubmission)), rep(0, NROW(farmsubmission)))
)

## ----estimatePopsizeFit-7-----------------------------------------------------
data.frame(IRLS  = round(c(res$beta, -ll(res$beta), res$iter), 4),
           optim = round(c(res2$beta, -ll(res2$beta), res2$iter[1]), 4))

## ----custom-family-example----------------------------------------------------
myFamilyFunction <- function(lambdaLink = c("logit", "cloglog", "probit"),
                             piLink     = c("logit", "cloglog", "probit"),
                             ...) {
  if (missing(lambdaLink)) lambdaLink <- "logit"
  if (missing(piLink))         piLink <- "logit"

  links <- list()
  attr(links, "linkNames") <- c(lambdaLink, piLink)

  lambdaLink <- switch(lambdaLink,
                       "logit"   = singleRcapture:::singleRinternallogitLink,
                       "cloglog" = singleRcapture:::singleRinternalcloglogLink,
                       "probit"  = singleRcapture:::singleRinternalprobitLink
  )

  piLink <- switch(piLink,
                   "logit"   = singleRcapture:::singleRinternallogitLink,
                   "cloglog" = singleRcapture:::singleRinternalcloglogLink,
                   "probit"  = singleRcapture:::singleRinternalprobitLink
  )

  links[1:2] <- c(lambdaLink, piLink)

  mu.eta <- function(eta, type = "trunc", deriv = FALSE, ...) {
    pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
    lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2

    if (!deriv) {
      switch (type,
              "nontrunc" = pi + 2 * lambda,
              "trunc" = 1 + lambda / (pi + lambda)
      )
    } else {
      # Only necessary if one wishes to use standard errors in predict method
      switch (type,
              "nontrunc" = {
                matrix(c(2, 1) * c(
                  lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / 2,
                  piLink(eta[, 2], inverse = TRUE, deriv = 1) / 2
                ), ncol = 2)
              },
              "trunc" = {
                matrix(c(
                  pi / (pi + lambda) ^ 2,
                  -lambda / (pi + lambda) ^ 2
                ) * c(
                  lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / 2,
                  piLink(eta[, 2], inverse = TRUE, deriv = 1) / 2
                ), ncol = 2)
              }
      )
    }
  }

  variance <- function(eta, type = "nontrunc", ...) {
    pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
    lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2

    switch (type,
            "nontrunc" = pi * (1 - pi) + 4 * lambda * (1 - lambda - pi),
            "trunc" = lambda * (1 - lambda) / (pi + lambda)
    )
  }

  Wfun <- function(prior, y, eta, ...) {
    pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
    lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2

    G01 <- ((lambda + pi) ^ (-2)) * piLink(eta[, 2], inverse = TRUE, deriv = 1) *
      lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) * prior / 4

    G00 <- ((lambda + pi) ^ (-2)) - (pi ^ (-2)) - lambda / ((lambda + pi) * (pi ^ 2))
    G00 <- G00 * prior * (piLink(eta[, 2], inverse = TRUE, deriv = 1) ^ 2) / 4

    G11 <- ((lambda + pi) ^ (-2)) - (((lambda + pi) * lambda) ^ -1)
    G11 <- G11 * prior * (lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) ^ 2) / 4

    matrix(
      -c(G11, # lambda
         G01, # mixed
         G01, # mixed
         G00  # pi
      ),
      dimnames = list(rownames(eta), c("lambda", "mixed", "mixed", "pi")),
      ncol = 4
    )
  }

  funcZ <- function(eta, weight, y, prior, ...) {
    pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
    lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2

    weight <- weight / prior

    G0 <- (2 - y) / pi     - ((lambda + pi) ^ -1)
    G1 <- (y - 1) / lambda - ((lambda + pi) ^ -1)

    G1 <- G1 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / 2
    G0 <- G0 *     piLink(eta[, 2], inverse = TRUE, deriv = 1) / 2

    uMatrix <- matrix(c(G1, G0), ncol = 2)

    weight <- lapply(X = 1:nrow(weight), FUN = function (x) {
      matrix(as.numeric(weight[x, ]), ncol = 2)
    })

    pseudoResid <- sapply(X = 1:length(weight), FUN = function (x) {
      #xx <- chol2inv(chol(weight[[x]])) # less computationally demanding
      xx <- solve(weight[[x]]) # more stable
      xx %*% uMatrix[x, ]
    })
    pseudoResid <- t(pseudoResid)
    dimnames(pseudoResid) <- dimnames(eta)
    pseudoResid
  }

  minusLogLike <- function(y, X, offset,
                           weight    = 1,
                           NbyK      = FALSE,
                           vectorDer = FALSE,
                           deriv     = 0,
                           ...) {
    y <- as.numeric(y)
    if (is.null(weight)) {
      weight <- 1
    }
    if (missing(offset)) {
      offset <- cbind(rep(0, NROW(X) / 2), rep(0, NROW(X) / 2))
    }

    if (!(deriv %in% c(0, 1, 2)))
      stop("Only score function and derivatives up to 2 are supported.")
    deriv <- deriv + 1

    switch (deriv,
            function(beta) {
              eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
              pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
              lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2
              -sum(weight * ((2 - y) * log(pi) + (y - 1) * log(lambda) - log(pi + lambda)))
            },
            function(beta) {
              eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
              pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
              lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2

              G0 <- (2 - y) / pi     - ((lambda + pi) ^ -1)
              G1 <- (y - 1) / lambda - ((lambda + pi) ^ -1)

              G1 <- G1 * weight * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / 2
              G0 <- G0 * weight *     piLink(eta[, 2], inverse = TRUE, deriv = 1) / 2

              if (NbyK) {
                XX <- 1:(attr(X, "hwm")[1])
                return(cbind(as.data.frame(X[1:nrow(eta), XX]) * G1,
                             as.data.frame(X[-(1:nrow(eta)), -XX]) * G0))
              }
              if (vectorDer) {
                return(cbind(G1, G0))
              }

              as.numeric(c(G1, G0) %*% X)
            },
            function (beta) {
              lambdaPredNumber <- attr(X, "hwm")[1]
              eta <- matrix(as.matrix(X) %*% beta, ncol = 2) + offset
              pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
              lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2

              res <- matrix(nrow = length(beta), ncol = length(beta),
                            dimnames = list(names(beta), names(beta)))

              # pi^2 derivative
              dpi <- (2 - y) / pi - (lambda + pi) ^ -1
              G00 <- ((lambda + pi) ^ (-2)) - (2 - y) / (pi ^ 2)

              G00 <- t(as.data.frame(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)] *
                                       (G00 * ((piLink(eta[, 2], inverse = TRUE, deriv = 1) / 2) ^ 2) +
                                          dpi * piLink(eta[, 2], inverse = TRUE, deriv = 2) / 2) * weight)) %*%
                as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
              # mixed derivative
              G01 <- (lambda + pi) ^ (-2)

              G01 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber]) *
                         G01 * (lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / 2) *
                         (piLink(eta[, 2], inverse = TRUE, deriv = 1) / 2) * weight) %*%
                as.matrix(X[-(1:(nrow(X) / 2)), -(1:lambdaPredNumber)])
              # lambda^2 derivative
              G11 <- ((lambda + pi) ^ (-2)) - (y - 1) / (lambda ^ 2)
              dlambda <- (y - 1) / lambda - ((lambda + pi) ^ -1)

              G11 <- t(as.data.frame(X[1:(nrow(X) / 2), 1:lambdaPredNumber] *
                                       (G11 * ((lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / 2) ^ 2) +
                                          dlambda * lambdaLink(eta[, 1], inverse = TRUE, deriv = 2) / 2) * weight)) %*%
                X[1:(nrow(X) / 2), 1:lambdaPredNumber]

              res[-(1:lambdaPredNumber), -(1:lambdaPredNumber)] <- G00
              res[1:lambdaPredNumber, 1:lambdaPredNumber] <- G11
              res[1:lambdaPredNumber, -(1:lambdaPredNumber)] <- t(G01)
              res[-(1:lambdaPredNumber), 1:lambdaPredNumber] <- G01

              res
            }
    )
  }

  validmu <- function(mu) {
    (sum(!is.finite(mu)) == 0) && all(0 < mu) && all(2 > mu)
  }

  # this is optional
  devResids <- function(y, eta, wt, ...) {
    0
  }

  pointEst <- function (pw, eta, contr = FALSE, ...) {
    pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
    lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2
    N <- pw / (lambda + pi)
    if(!contr) {
      N <- sum(N)
    }
    N
  }

  popVar <- function (pw, eta, cov, Xvlm, ...) {
    pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
    lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2

    bigTheta1 <- -pw / (pi + lambda) ^ 2 # w.r to pi
    bigTheta1 <- bigTheta1 * piLink(eta[, 2], inverse = TRUE, deriv = 1) / 2
    bigTheta2 <- -pw / (pi + lambda) ^ 2 # w.r to lambda
    bigTheta2 <- bigTheta2 * lambdaLink(eta[, 1], inverse = TRUE, deriv = 1) / 2 # w.r to lambda

    bigTheta <- t(c(bigTheta2, bigTheta1) %*% Xvlm)

    f1 <- t(bigTheta) %*% as.matrix(cov) %*% bigTheta

    f2 <- sum(pw * (1 - pi - lambda) / ((pi + lambda) ^ 2))

    f1 + f2
  }

  dFun <- function (x, eta, type = c("trunc", "nontrunc")) {
    if (missing(type)) type <- "trunc"
    pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
    lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2

    switch (type,
            "trunc" = {
              (pi * as.numeric(x == 1) + lambda * as.numeric(x == 2)) / (pi + lambda)
            },
            "nontrunc" = {
              (1 - pi - lambda) * as.numeric(x == 0) +
                pi * as.numeric(x == 1) + lambda * as.numeric(x == 2)
            }
    )
  }

  simulate <- function(n, eta, lower = 0, upper = Inf) {
    pi     <-     piLink(eta[, 2], inverse = TRUE) / 2
    lambda <- lambdaLink(eta[, 1], inverse = TRUE) / 2
    CDF <- function(x) {
      ifelse(x == Inf, 1,
             ifelse(x < 0, 0,
                    ifelse(x < 1, 1 - pi - lambda,
                           ifelse(x < 2, 1 - lambda, 1))))
    }
    lb <- CDF(lower)
    ub <- CDF(upper)
    p_u <- stats::runif(n, lb, ub)
    sims <- rep(0, n)
    cond <- CDF(sims) <= p_u
    while (any(cond)) {
      sims[cond] <- sims[cond] + 1
      cond <- CDF(sims) <= p_u
    }
    sims
  }

  getStart <- expression(
    if (method == "IRLS") {
      etaStart <- cbind(
        family$links[[1]](mean(observed == 2) * (1 + 0 * (observed == 2))), # lambda
        family$links[[2]](mean(observed == 1) * (1 + 0 * (observed == 1)))  # pi
      ) + offset
    } else if (method == "optim") {
      init <- c(
        family$links[[1]](weighted.mean(observed == 2, priorWeights) * 1 + .0001),
        family$links[[2]](weighted.mean(observed == 1, priorWeights) * 1 + .0001)
      )
      if (attr(terms, "intercept")) {
        coefStart <- c(init[1], rep(0, attr(Xvlm, "hwm")[1] - 1))
      } else {
        coefStart <- rep(init[1] / attr(Xvlm, "hwm")[1], attr(Xvlm, "hwm")[1])
      }
      if ("(Intercept):pi" %in% colnames(Xvlm)) {
        coefStart <- c(coefStart, init[2], rep(0, attr(Xvlm, "hwm")[2] - 1))
      } else {
        coefStart <- c(coefStart, rep(init[2] / attr(Xvlm, "hwm")[2], attr(Xvlm, "hwm")[2]))
      }
    }
  )

  structure(
    list(
      makeMinusLogLike = minusLogLike,
      densityFunction  = dFun,
      links     = links,
      mu.eta    = mu.eta,
      valideta  = function (eta) {TRUE},
      variance  = variance,
      Wfun      = Wfun,
      funcZ     = funcZ,
      devResids = devResids,
      validmu   = validmu,
      pointEst  = pointEst,
      popVar    = popVar,
      family    = "myFamilyFunction",
      etaNames  = c("lambda", "pi"),
      simulate  = simulate,
      getStart  = getStart,
      extraInfo = c(
        mean       = "pi / 2 + lambda",
        variance   = paste0("(pi / 2) * (1 - pi / 2) + 2 * lambda * (1 - lambda / 2 - pi / 2)"),
        popSizeEst = "(1 - (pi + lambda) / 2) ^ -1",
        meanTr     = "1 + lambda / (pi + lambda)",
        varianceTr = paste0("lambda * (1 - lambda / 2) / (pi + lambda)")
      )
    ),
    class = c("singleRfamily", "family")
  )
}

## -----------------------------------------------------------------------------
set.seed(123)
Y <- simulate(
  myFamilyFunction(lambdaLink = "logit", piLink = "logit"),
  nsim = 1000, eta = matrix(0, nrow = 1000, ncol = 2),
  truncated = FALSE
)
mm <- estimatePopsize(
  formula = Y ~ 1,
  data = data.frame(Y = Y[Y > 0]),
  model = myFamilyFunction(lambdaLink = "logit",
                           piLink = "logit"),
  # the usual observed information matrix
  # is ill-suited for this distribution
  controlPopVar = controlPopVar(covType = "Fisher")
)
summary(mm)

## ----single-link--------------------------------------------------------------
singleRcapture:::singleRinternalcloglogLink

