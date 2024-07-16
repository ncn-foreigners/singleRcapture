#N <- 1000
###gender <- rbinom(N, 1, 0.2)
gender <- rep(0:1, c(812, 188))

etaLambda <- -1   + 0.5 * gender
etaAlpha  <- -1.7 + 1 * gender
etaOmega  <- -.2  + 0.8 * gender
etaPi     <- -.3  + 0.9 * gender

expect_silent(
  res <- data.frame(
    ztpoisson = simulate(
      ztpoisson(), 
      eta = cbind(etaLambda), 
      seed = 1
    ),
    ztgeom = simulate(
      ztgeom(), 
      eta = cbind(etaLambda), 
      seed = 1
    ),
    ztnegbin = simulate(
      ztnegbin(), 
      eta = cbind(etaLambda, etaAlpha), 
      seed = 1
    ),
    zotpoisson = simulate(
      zotpoisson(), 
      eta = cbind(etaLambda), 
      seed = 1
    ),
    zotgeom = simulate(
      zotgeom(), 
      eta = cbind(etaLambda), 
      seed = 1
    ),
    zotnegbin = simulate(
      zotnegbin(), 
      eta = cbind(etaLambda, etaAlpha), 
      seed = 1
    ),
    ztoipoisson = simulate(
      ztoipoisson(), 
      eta = cbind(etaLambda, etaOmega), 
      seed = 1
    ),
    ztoigeom = simulate(
      ztoigeom(), 
      eta = cbind(etaLambda, etaOmega), 
      seed = 1
    ),
    ztoinegbin = simulate(
      ztoinegbin(), 
      eta = cbind(etaLambda, etaAlpha, etaOmega), 
      seed = 1
    ),
    oiztpoisson = simulate(
      oiztpoisson(), 
      eta = cbind(etaLambda, etaOmega), 
      seed = 1
    ),
    oiztgeom = simulate(
      oiztgeom(), 
      eta = cbind(etaLambda, etaOmega), 
      seed = 1
    ),
    oiztnegbin = simulate(
      oiztnegbin(), 
      eta = cbind(etaLambda, etaAlpha, etaOmega), 
      seed = 1
    ),
    ztHurdlepoisson = simulate(
      ztHurdlepoisson(), 
      eta = cbind(etaLambda, etaPi), 
      seed = 1
    ),
    ztHurdlegeom = simulate(
      ztHurdlegeom(), 
      eta = cbind(etaLambda, etaPi), 
      seed = 1
    ),
    ztHurdlenegbin = simulate(
      ztHurdlenegbin(), 
      eta = cbind(etaLambda, etaAlpha, etaPi), 
      seed = 1
    ),
    Hurdleztpoisson = simulate(
      Hurdleztpoisson(), 
      eta = cbind(etaLambda, etaPi), 
      seed = 1
    ),
    Hurdleztgeom = simulate(
      Hurdleztgeom(), 
      eta = cbind(etaLambda, etaPi), 
      seed = 1
    ),
    Hurdleztnegbin = simulate(
      Hurdleztnegbin() , 
      eta = cbind(etaLambda, etaAlpha, etaPi), 
      seed = 1
    ),
    chao = simulate(
      chao(), 
      eta = cbind(etaLambda), 
      seed = 1
    ),
    zelterman = simulate(
      zelterman(), 
      eta = cbind(etaLambda), 
      seed = 1
    )
  )
)

expect_true(
  all(res$chao == res$zelterman)
)

expect_true(
  all(res$zotpoisson == res$ztpoisson)
)

expect_true(
  all(res$zotgeom == res$ztgeom)
)

expect_true(
  all(res$zotnegbin == res$ztnegbin)
)

