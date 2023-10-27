#-------------------------------------------------------------------------------
#                 Oswaldo Gressani, Metropolis example.
#                 Copyright, 2023. All rights reserved.
#------------------------------------------------------------------------------

# Package to sample from multivariate normal
library("mvtnorm")
library("coda")


Metropolis <- function(M = 31000, burn = 1000, seed = 123,
                       alpha_start=-36.95, beta_start = 20.87){

# Seed
set.seed(seed)
    
# Encode data
w_i <- c(1.583, 1.712, 1.774, 1.843, 1.875, 1.892,1.902, 1.930)
y_i <- c(7, 12, 18, 50, 59, 60, 61, 64)
n_i <- c(58, 61, 63, 55, 61, 68, 63, 64)

# Wise starting values (frequentist approach)
freqfit <- glm(cbind(y_i, n_i-y_i)~w_i, family = "binomial")
Sigma_hat <- vcov(freqfit)
alpha_hat <- alpha_start
beta_hat <- beta_start


# Joint log posterior distribution
lpost <- function(theta) {
  eta <- theta[1] + theta[2] * w_i
  pi <- exp(eta) / (1 + exp(eta))
  return(sum(dbinom(y_i, n_i, pi, log = TRUE)))
}

# Starting value for alpha and beta
M <- 31000
burnin <- burn
theta <- array(dim = c(2, M))
theta[, 1] <- c(alpha_hat, beta_hat)
n_accept <- 0
sd_prop <- 2.4^2

# Metropolis loop
for (i in 2:M){
  theta_prop <- theta[, i - 1] + rmvnorm(n = 1, mean = c(0, 0), 
                                         sigma = sd_prop * Sigma_hat)
  logr <- lpost(theta_prop) - lpost(theta[,(i-1)])
  u <- runif(1)
  if(logr >= 0 || u <= exp(logr)){
    theta[,i] <- theta_prop
    n_accept <- n_accept + 1
  } else {
    theta[,i] <- theta[,i-1]
  }
}

# Exclude burnin
theta <- theta[,-c(1:burnin)]
accept_rate <- round(n_accept/(M-1),digits=2) * 100
rownames(theta) <- c("alpha", "beta")

outlist <- list(theta = theta,
                M=M,
                burnin = burnin,
                accept_rate = accept_rate,
                lpost = lpost
                )

return(outlist)

}

runmean <- function(x){
  nlen <- length(x)
  mean_iteration <- c()
  mean_iteration[1] <- x[1]
  for(j in 2:nlen){
    mean_iteration[j] <- mean(x[1:j])
  }
  return(mean_iteration)
}

#------------------------- Assessing convergence -------------------------------

MCMCrun <- Metropolis()
lpost <- MCMCrun$lpost
theta <- MCMCrun$theta
M <- MCMCrun$M
burnin <- MCMCrun$burnin

#--- Basic summary 
summary(as.mcmc(theta[1,]))

#--- Trace plots of alpha and beta
par(mfrow = c(1,2))
coda::traceplot(as.mcmc(theta[1,]), main = "alpha", col = "firebrick1")
coda::traceplot(as.mcmc(theta[2,]), main = "beta", col = "darkgreen")

#--- Trace plot of the logposterior (lpost)
lpost_chain <- as.mcmc(apply(t(theta), MARGIN = 1, FUN = lpost))
par(mfrow = c(1,1))
coda::traceplot(lpost_chain, ylab="Log(posterior)")


#--- Autocorrelation plots
coda::autocorr.plot(as.mcmc(theta[1,]), lag.max = 30, main = "alpha", lwd = 2)
coda::autocorr(as.mcmc(theta[1,]), lags = c(2,10)) # Specific lags
coda::autocorr.plot(as.mcmc(theta[2,]), lag.max = 30, main = "beta", lwd = 2)
coda::autocorr.plot(lpost_chain, lag.max = 30, main = "LOGPOST", lwd = 2)


#--- Ergodic mean plots
ergtheta1 <- runmean(theta[1,])
ergtheta2 <- runmean(theta[2,])
par(mfrow = c(2,1))
plot(seq_len(M-burnin), ergtheta1, type = "l", xlab = "Iterations", 
     ylab = "alpha", col = "blue", lwd = 2)
plot(seq_len(M-burnin), ergtheta2, type = "l", xlab = "Iterations", 
     ylab = "beta", col = "firebrick", lwd = 2)

#--- QQ-plots (first half and second half of alpha)
par(mfrow = c(1,1))
alpha_1st <- theta[1,(1: (M-burnin)/2)]
alpha_2nd <- theta[1,(((M-burnin)/2 + 1):(M-burnin))]
qqplot(alpha_1st, alpha_2nd, xlab = "First half of alpha",
       ylab = "Second half of alpha", col = "red")

#--- Geweke diagnostic
coda::geweke.diag(as.mcmc(theta[1,]))
coda::geweke.diag(as.mcmc(theta[2,]))
# head(as.mcmc(t(theta)))
coda::geweke.plot(as.mcmc(t(theta)))

#--- Brooks-Gelman-Rubin (BGR) diagnostic

#-- Generate 3 chains for each variable
MCMCrun1 <- Metropolis(seed=123,  alpha_start=-36.95, beta_start = 20.87)
MCMCrun2 <- Metropolis(seed=1990, alpha_start=-45, beta_start = 15)
MCMCrun3 <- Metropolis(seed=2020, alpha_start=-25, beta_start = 12)

MCMCrun1$accept_rate
MCMCrun2$accept_rate
MCMCrun3$accept_rate

theta1 <- MCMCrun1$theta
theta2 <- MCMCrun2$theta
theta3 <- MCMCrun3$theta

#--- Visualize trace plots for alpha and beta
par(mfrow = c(1,2))
#-- alpha
plot(seq_len(M-burnin), theta1[1,], type = "l", col = "firebrick1",
     xlab = "Iterations", ylab = "", main = "alpha")
lines(seq_len(M-burnin), theta2[1,], type = "l", col = "green")
lines(seq_len(M-burnin), theta3[1,], type = "l", col = "darkorange1")
#-- beta
plot(seq_len(M-burnin), theta1[2,], type = "l", col = "black",
     xlab = "Iterations", ylab = "", main = "beta")
lines(seq_len(M-burnin), theta2[2,], type = "l", col = "dodgerblue2")
lines(seq_len(M-burnin), theta3[2,], type = "l", col = "gold1")

#--- BGR diagnostic
alpha.obj <- mcmc.list(chain1 = mcmc(theta1[1,]), chain2 = mcmc(theta2[1,]),
                       chain3 = mcmc(theta3[1,]))
beta.obj <- mcmc.list(chain1 = mcmc(theta1[2,]), chain2 = mcmc(theta2[2,]),
                       chain3 = mcmc(theta3[2,]))
coda::gelman.plot(alpha.obj, main = "alpha")
coda::gelman.plot(beta.obj, main = "beta")


