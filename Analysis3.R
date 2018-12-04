
library(sandwich)
data(PublicSchools); ps <- na.omit(PublicSchools)
ps$Income <- ps$Income * 0.01
#Data 
ps[, 'Income2'] <- ps[,'Income']^2; ps[, 'Intercept'] <- 1
ps <- ps[c("Expenditure", "Intercept", "Income", "Income2")]

# Y matrix
Y. <- as.matrix.data.frame(ps[c(1)]); Y. <- Y.[-2, ]
N <- dim(X.)[1]


# Standard deviation of Y for currupting data 
original <- sd(Y.)

#Iterations
K <- 5000

#Store distribution
LAMBDA <- BETA <- NULL
S2 <- rep(NA, K)

#Prior Values
v <- 4 #df 

#Current value placeholders
lambda <- rep(NA, N)
RSS <- rep(NA, N)

library(MASS); library(invgamma); library(MCMCpack); library(mvtnorm)

#Robust Regression
RobustReg <- function(X, Y, beta.0. , sigma2.0. , N. = N, K. = K, v. = v, RSS. = RSS, lambda. = lambda, LAMBDA. = LAMBDA, S2. = S2, BETA. = BETA ) {
  for (i in 1:K.) {
    if (i == 1){
      beta <- beta.0.
      sigma2 <- sigma2.0.
    }
    #lambdas
    for (n in 1:N.){
      rss.i <- (Y[n] - t(beta) %*% X[n, ])^2
      RSS.[n] <- rss.i
      lambda.[n] <- rinvgamma(1, (v. +1)/2, (v. + ((sigma2)^(-1))*rss.i)/2)
    }
    lambda.matrix <- diag(lambda., nrow = N., ncol = N.)
    LAMBDA. <- rbind(lambda., LAMBDA.)
    
    #sigma2
    wRSS <- sum(RSS./lambda.)
    sigma2 <- rinvgamma(1, (N. + 2*(1 - 1))/2, (1/2)*wRSS)
    S2.[i] <- sigma2
    
    #Beta 
    beta.sigma <- solve(t(X)%*% solve(sigma2*lambda.matrix) %*% X )
    beta.mu <- beta.sigma %*% (t(X) %*% solve(sigma2*lambda.matrix) %*%Y  )
    beta <- as.matrix(mvrnorm(1, beta.mu, beta.sigma))
    BETA. <- rbind(t(beta),BETA.)
  }
  
  return(list(BETA = BETA., S2 = S2., LAMBDA = LAMBDA.))
}


Corrupt2 <- function(p, dt,  X_ , Y_ = Y.,  beta.0_ , sigma2.0_ , 
                     N. = N, K. = K, v. = v, RSS. = RSS, lambda. = lambda, LAMBDA. = LAMBDA, 
                     S2. = S2, BETA. = BETA)
{
  #corruption 
  corrupt <- which(dt[,"Income"] > 89)
  noise <- rexp(sum(with(dt, Income > 89)), p) # generate the noise to add
  if (p == 0){ noise = 0 }
  Y_[corrupt] <- Y_[corrupt] + noise  
  reg <- RobustReg(X = X_, Y = Y_, beta.0. = beta.0_, sigma2.0. = sigma2.0_)
  
  #Gaussian
  J <- dim(X_)[2]
  dt[corrupt,"Expenditure"] <- dt[corrupt, "Expenditure"] + noise
  m <- lm(Expenditure ~ Income, data = dt)
  if (J == 3){ m <- lm(Expenditure ~ Income + I(Income^2), data = dt) }
  ols <- as.array(m$coefficients)
  df <- N. - 2
  s2 <- sum((m$residuals)^2)
  l.BETA.st <- rmvt(K., sigma = s2*solve((t(X_) %*% X_)), df =df, delta = ols )
  l.S2.st <- rinvgamma(K., df/2, s2/2)
  
  b <- as.array(apply(reg$BETA, MARGIN = 2, FUN = mean))
  gb <- as.array(apply(l.BETA.st, MARGIN = 2, FUN = mean))
  mse <- mean((Y_ - X_ %*% b)^2)
  gmse <- mean((Y_ - X_ %*% gb)^2)
  
  return(list(BETA = reg$BETA, S2 = reg$S2, LAMBDA = reg$LAMBDA, gBETA = l.BETA.st, gS2 = l.S2.st, mse = mse, 
              gmse = gmse, rb = b, gb = gb))
}

#Linear, drop alaska
lX. <- as.matrix.data.frame(ps[c("Intercept", "Income")]); lX. <- lX.[-2,]
#quadratic 
X. <- as.matrix(ps[c("Intercept", "Income", "Income2")]); X. <- X.[-2, ]


# Starting values
#Linear
ld <- ps[c("Expenditure", "Intercept", "Income")]; ld <- d[-2,]
m2 <- lm(Expenditure ~ Income , data = ld)
lbeta.0 <- as.matrix(m2$coefficients)
lsigma2.0 <- summary(m2)$sigma

#quadratic 
qd <- ps[-2, ]
m3 <- lm(Expenditure ~ Income + I(Income^2), data = qd)
beta.0 <- as.matrix(m3$coefficients)
sigma2.0 <- summary(m3)$sigma

exp5 <- Corrupt2(p = 0.001, dt = qd, X_ = X. , beta.0_ = beta.0, sigma2.0_ = sigma2.0)
exp5_l <- Corrupt2(p = 0.001, dt = ld, X_ = lX. , beta.0_ = lbeta.0, sigma2.0_ = lsigma2.0)
exp4 <- Corrupt2(p = 0.0025, dt = qd, X_ = X. , beta.0_ = beta.0, sigma2.0_ = sigma2.0)
exp4_l <- Corrupt2(p = 0.0025, dt = ld, X_ = lX. , beta.0_ = lbeta.0, sigma2.0_ = lsigma2.0)
exp3 <- Corrupt2(p = 0.005, dt = qd, X_ = X. , beta.0_ = beta.0, sigma2.0_ = sigma2.0)
exp3_l <- Corrupt2(p = 0.005, dt = ld, X_ = lX. , beta.0_ = lbeta.0, sigma2.0_ = lsigma2.0)
exp2 <- Corrupt2(p = 0.01, dt = qd, X_ = X. , beta.0_ = beta.0, sigma2.0_ = sigma2.0)
exp2_l <- Corrupt2(p = 0.01, dt = ld, X_ = lX. , beta.0_ = lbeta.0, sigma2.0_ = lsigma2.0)
exp1 <- Corrupt2(p = 0.03, dt = qd, X_ = X. , beta.0_ = beta.0, sigma2.0_ = sigma2.0)
exp1_l <- Corrupt2(p = 0.03, dt = ld, X_ = lX. , beta.0_ = lbeta.0, sigma2.0_ = lsigma2.0)


save(exp1, exp1_l, exp2, exp2_l, exp3, exp3_l, exp4, exp4_l, exp5, exp5_l, file="/Users/MacUser/Desktop/MA578/Final Project/Corrupt2alt.RData")


