# test simulate
set.seed(123)
expect_equivalent(
  {N <- 10000
   gender <- rbinom(N, 1, 0.2)
   eta <- -1 + 0.5*gender
   disp <- 1
   counts <- rnbinom(N, mu = exp(eta), size = disp)
   df <- data.frame(gender, eta, counts)
   df2 <- subset(df, counts > 0)
   mod1 <-  estimate_popsize(formula = counts ~ 1 + gender, 
                             data = df2,
                             model = "ztnegbin", 
                             method = "mle",
                             pop.var = "analytic")
   mid1_sim <- simulate(mod1, 10)
   dim(mid1_sim)
  },
  c(2920, 10)
)