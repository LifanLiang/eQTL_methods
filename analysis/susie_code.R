library(susieR)

X <- t(readRDS("data/neuron_stim_genotype_chr9_2e6-3e6.rds")[,-1])

pi <- rep(c(2,4,8,16,32) * 1e-4, 673)

sim_group_pi <- function(X, pi=pi, sigma=0.25, sigma.e=1.0) {
  ### Impose the constraint that there is at least one causal SNP
  N <- 0
  while(N==0) {
    gamma <- rbinom(ncol(X),1,pi)
    N <- sum(gamma)
  }
  beta <- rnorm(ncol(X),sd=sigma) * gamma
  eps <- rnorm(nrow(X),sd=sigma.e)
  X <- scale(X, center=T, scale=T)
  Y <- drop(X %*% beta + eps)
  list(y=Y, X=X, pi=pi, gamma=gamma, beta=beta, eps=eps)
}

set.seed(22)
sim.res <- sim_group_pi(X,pi)

susie.fit <- susie(sim.res$X, sim.res$y, prior_weights = sim.res$pi,
                   L=sum(sim.res$gamma), coverage = 0.95, standardize = F,
                   intercept = T, compute_univariate_zscore = T)

plot(predict(susie.fit), sim.res$y)