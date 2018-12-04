
library(sandwich)
data(PublicSchools); ps <- na.omit(PublicSchools)
ps$Income <- ps$Income * 0.01
#Data 
ps[, 'Income2'] <- ps[,'Income']^2; ps[, 'Intercept'] <- 1
ps <- ps[c("Expenditure", "Intercept", "Income", "Income2")]

#Linear, drop alaska
X. <- as.matrix.data.frame(ps[c("Intercept", "Income")]); X. <- X.[-2,]
Y. <- as.matrix.data.frame(ps[c(1)]); Y. <- Y.[-2, ]
N <- dim(X.)[1]

# Standard deviation of Y for currupting data 
original <- sd(Y.)

library(MASS); library(invgamma); library(MCMCpack)

#Iterations
K <- 5000

#Store distribution
LAMBDA <- BETA <- NULL
S2 <- rep(NA, K)

#Prior Values
v <- 4 #df 

# Starting values 
d <- ps[c("Expenditure", "Intercept", "Income")]; d <- d[-2,]
m2 <- lm(Expenditure ~ Income , data = d)
beta.0 <- as.matrix(m2$coefficients)
sigma2.0 <- summary(m2)$sigma


#Current value placeholders
lambda <- rep(NA, N)
RSS <- rep(NA, N)

#Robust Regression
set.seed(143)
RobustReg <- function(X, Y, beta.0. = beta.0, sigma2.0. = sigma2.0, N. = N, K. = K, v. = v, RSS. = RSS, lambda. = lambda, LAMBDA. = LAMBDA, S2. = S2, BETA. = BETA ) {
  for (i in 1:K.) {
    if (i == 1){
      beta <- beta.0.
      sigma2 <- sigma2.0
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


set.seed(2)
Corrupt <- function(p, dt = d, sd = original, X_ = X., Y_=Y.,  beta.0. = beta.0, sigma2.0. = sigma2.0, 
                    N. = N, K. = K, v. = v, RSS. = RSS, lambda. = lambda, LAMBDA. = LAMBDA, 
                    S2. = S2, BETA. = BETA)
  {
  #corruption 
  corrupt <- rbinom(N.,1,0.10)    # choose an average of 10% to corrupt at random
  corrupt <- as.logical(corrupt)
  noise <- rnorm(sum(corrupt), 0, sqrt(p)*sd) # generate the noise to add
  if (p == 0){ noise = 0 }
  Y_[corrupt] <- Y_[corrupt] + noise  
  reg <- RobustReg(X = X_, Y = Y_)
  
  #Gaussian
  dt[corrupt,"Expenditure"] <- dt[corrupt, "Expenditure"] + noise
  m <- lm(Expenditure ~ Income, data = dt)
  ols <- as.array(m$coefficients)
  df <- N. - 2
  s2 <- sum((m$residuals)^2)
  l.BETA.st <- rmvt(K., sigma = s2*solve((t(X_) %*% X_)), df =df, delta = ols )
  l.S2.st <- rinvgamma(K., df/2, s2/2)
  return(list(BETA = reg$BETA, S2 = reg$S2, LAMBDA = reg$LAMBDA, gBETA = l.BETA.st, gS2 = l.S2.st))
}

# Standard Regression
#Linear 
Gaussian <- function(model, N. = N, K. = K, X = X., Y = Y.) {
  ols.l.beta <- as.array(model$coefficients)
  df <- N. - 2 
  s2 <- sum((model$residuals)^2)
  l.BETA.st <- rmvt(K., sigma = s2*solve((t(X) %*% X)), df =df, delta = ols.l.beta )
  l.S2.st <- rinvgamma(K., df/2, s2/2)
  return (list(gBETA = l.BETA.st, gS2 = l.S2.st))
}

r0 <- Corrupt(p = 0 )
r1 <- Corrupt(p = 0.25)
r2 <- Corrupt(p = 0.5)
r3 <- Corrupt(p = 0.75)
r4 <- Corrupt(p = 1)
r5 <- Corrupt(p = 2)
r6 <- Corrupt(p = 3)
r7 <- Corrupt(p = 4)


Corrupt2 <- function(p, dt = d, sd = original, X_ = X., Y_=Y.,  beta.0. = beta.0, sigma2.0. = sigma2.0, 
                    N. = N, K. = K, v. = v, RSS. = RSS, lambda. = lambda, LAMBDA. = LAMBDA, 
                    S2. = S2, BETA. = BETA)
{
  #corruption 
  corrupt <- which(dt[,"Income"] > 89)
  noise <- rexp(sum(with(dt, Income > 89)), p) # generate the noise to add
  if (p == 0){ noise = 0 }
  Y_[corrupt] <- Y_[corrupt] + noise  
  reg <- RobustReg(X = X_, Y = Y_)
  
  #Gaussian
  dt[corrupt,"Expenditure"] <- dt[corrupt, "Expenditure"] + noise
  m <- lm(Expenditure ~ Income, data = dt)
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
              gmse = gmse))
}

save(r0, r1, r2, r3, r4, r5, r6, r7, file="/Users/MacUser/Desktop/MA578/Final Project/Corrupt.RData")


