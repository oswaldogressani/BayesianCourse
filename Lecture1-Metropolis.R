#-------------------------------------------------------------------------------
#                 Oswaldo Gressani, Metropolis example.
#                 Copyright, 2023. All rights reserved.
#                 Code inspired from Philippe Lambert, Liège University.
#-------------------------------------------------------------------------------

# Package to sample from multivariate normal
library("mvtnorm")
library("coda")

set.seed(123)

# Encode data
w_i <- c(1.583, 1.712, 1.774, 1.843, 1.875, 1.892,1.902, 1.930)
y_i <- c(7, 12, 18, 50, 59, 60, 61, 64)
n_i <- c(58, 61, 63, 55, 61, 68, 63, 64)

# Wise starting values (frequentist approach)
freqfit <- glm(cbind(y_i, n_i-y_i)~w_i, family = "binomial")
alpha_hat <- coef(freqfit)[1]
beta_hat <- coef(freqfit)[2]
Sigma_hat <- vcov(freqfit)

# Joint log posterior distribution
lpost <- function(theta) {
  eta <- theta[1] + theta[2] * w_i
  pi <- exp(eta) / (1 + exp(eta))
  return(sum(dbinom(y_i, n_i, pi, log = TRUE)))
}

# Starting value for alpha and beta
M <- 31000
burnin <- 1000
theta <- array(dim = c(2, M))
theta[, 1] <- c(alpha_hat, beta_hat)
n_accept <- 0
sd_prop <- 2.4

# Metropolis loop
for (i in 2:M){
  theta_prop <- theta[, i - 1] + rmvnorm(n = 1, mean = c(0, 0), 
                                         sigma = sd_prop ^ 2 * Sigma_hat)
  logr <- lpost(theta_prop) - lpost(theta[,(i-1)])
  u <- runif(1)
  if(logr >= 0 || u <= exp(logr)){
      theta[,i] <- theta_prop
      n_accept <- n_accept + 1
  } else {
    theta[,i] <- theta[,i-1]
  }
}

# Exclude warm up
theta <- theta[,-c(1:burnin)]
accept_rate <- round(n_accept/(M-1),digits=2) * 100
rownames(theta) <- c("alpha", "beta")

# Trace plots
par(mfrow = c(2,1))
coda::traceplot(as.mcmc(theta[1,]), main = "alpha", col = "orange")
coda::traceplot(as.mcmc(theta[2,]), main = "beta", col = "darkgreen")

par(mfrow = c(1,2))
#Histogram of generated samples for alpha
hist(theta[1,],main="alpha", xlab = "", freq = FALSE, ylim = c(0,0.13))
lines(density(theta[1,]), type = "l", col = "orange", lwd = 2)
HPDalpha <- HPDinterval(as.mcmc(t(theta)))[1, ]
CIalpha <- quantile(theta[1, ], probs = c(0.025, 0.975))
HPDalpha
CIalpha
abline(v = c(HPDalpha[1], HPDalpha[2]), col = "red", lty = 2)

#Histogram of generated samples for beta
hist(theta[2,],main="beta", xlab = "", freq = FALSE, ylim = c(0,0.25))
lines(density(theta[2,]), type = "l", col = "darkgreen", lwd = 2)
HPDbeta <- HPDinterval(as.mcmc(t(theta)))[2, ]
CIbeta <- quantile(theta[2, ], probs = c(0.025, 0.975))
HPDbeta
CIbeta
abline(v = c(HPDbeta[1], HPDbeta[2]), col = "red", lty = 2)

# Descriptive statistics for each model parameter
summary(t(theta))
cat("Acceptance rate --> ",accept_rate,"%.", "\n")

# Plot estimated probability of death
alpha_MCMC <- mean(theta[1,])
beta_MCMC <- mean(theta[2,])

par(mfrow = c(1,1))
wdom <- seq(1.5,2, length = 200)
pihat <- exp(alpha_hat + beta_hat * wdom) / (1 +exp(alpha_hat + beta_hat * wdom))
pihatMetro <- exp(alpha_MCMC + beta_MCMC * wdom) / 
  (1 + exp(alpha_MCMC + beta_MCMC * wdom))
plot(wdom, pihat, ylab = "Death prob of mice", xlab = "w (virus conc.)", 
     type = "l", lwd = 2, col = "green")
lines(wdom, pihatMetro, col = "red", lwd = 1, lty = 2)
grid(nx = 10, ny = 10)
legend("topleft", c("Frequentist","MCMC"), col = c("green","red"), lty = c(1,2),
       bty = "n")


