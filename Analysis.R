library(sandwich)
data(PublicSchools); ps <- na.omit(PublicSchools)
ps$Income <- ps$Income * 0.01

#Robust Regression
library(MASS); library(invgamma); library(MCMCpack)

K <- 5000
#Store distribution
LAMBDA <- BETA <- NULL
S2 <- rep(NA, K)

#Prior Values
v <- 4 #df 

#Jeffrey's prior
a <- c(1, 3/2)

#Quadratic
#Data 
ps[, 'Income2'] <- ps[,'Income']^2; ps[, 'Intercept'] <- 1
ps <- ps[c("Expenditure", "Intercept", "Income", "Income2")]
X <- as.matrix.data.frame(ps[-c(1)])
Y <- as.matrix.data.frame(ps[c(1)])
N <- dim(ps)[1]

# Starting value of BETA and S2
m1 <- lm(Expenditure ~ Income + I(Income^2), data = ps)
beta <- as.matrix(m1$coefficients)
sigma2 <- summary(m1)$sigma

#Current value 
lambda <- rep(NA, N)
RSS <- rep(NA, N)


#Robust Regression
set.seed(143)
for (i in 1:K) {
  #lambdas
  for (n in 1:N){
    rss.i <- (Y[n] - t(beta) %*% X[n, ])^2
    RSS[n] <- rss.i
    lambda[n] <- rinvgamma(1, (v +1)/2, (v + ((sigma2)^(-1))*rss.i)/2)
  }
  lambda.matrix <- diag(lambda, nrow = N, ncol = N)
  LAMBDA <- rbind(lambda, LAMBDA)
  
  #sigma2
  wRSS <- sum(RSS/lambda)
  sigma2 <- rinvgamma(1, (N + 2*(a[1] - 1))/2, (1/2)*wRSS)
  S2[i] <- sigma2
  
  #Beta 
  beta.sigma <- solve(t(X)%*% solve(sigma2*lambda.matrix) %*% X )
  beta.mu <- beta.sigma %*% (t(X) %*% solve(sigma2*lambda.matrix) %*%Y  )
  beta <- as.matrix(mvrnorm(1, beta.mu, beta.sigma))
  BETA <- rbind(t(beta),BETA)
}

#Standard Regression 
#Quadratic
ols.beta <- as.array(m1$coefficients)
library(mvtnorm)
df <- N + a[1] - 3 - 1
s2 <- sum((m1$residuals)^2)
BETA.st <- rmvt(K, sigma = s2*solve((t(X) %*% X)), df =df, delta = ols.beta )
S2.st <- rinvgamma(K, df/2, s2/2)

#With no quardratic term 

X <- as.matrix.data.frame(ps[-c(1, 4)])
Y <- as.matrix.data.frame(ps[c(1)])
N <- dim(ps)[1]

# Starting value of BETA and S2
m2 <- lm(Expenditure ~ Income, data = ps[-c(4)])
beta <- as.matrix(m2$coefficients)
sigma2 <- summary(m2)$sigma

#Store distribution
l.LAMBDA <- l.BETA <- NULL
l.S2 <- rep(NA, K)

set.seed(143)
for (i in 1:K) {
  #lambdas
  for (n in 1:N){
    rss.i <- (Y[n] - t(beta) %*% X[n, ])^2
    RSS[n] <- rss.i
    lambda[n] <- rinvgamma(1, (v +1)/2, (v + ((sigma2)^(-1))*rss.i)/2)
  }
  lambda.matrix <- diag(lambda, nrow = N, ncol = N)
  l.LAMBDA <- rbind(lambda, l.LAMBDA)
  
  #sigma2
  wRSS <- sum(RSS/lambda)
  sigma2 <- rinvgamma(1, (N + 2*(a[1] - 1))/2, (1/2)*wRSS)
  l.S2[i] <- sigma2
  
  #Beta 
  beta.sigma <- solve(t(X)%*% solve(sigma2*lambda.matrix) %*% X )
  beta.mu <- beta.sigma %*% (t(X) %*% solve(sigma2*lambda.matrix) %*%Y  )
  beta <- as.matrix(mvrnorm(1, beta.mu, beta.sigma))
  l.BETA <- rbind(t(beta),l.BETA)
}
library("coda")
effectiveSize(BETA)


# Standard Regression
#Linear 
ols.l.beta <- as.array(m2$coefficients)
df <- N + a[1] - 2 - 1
s2 <- sum((m2$residuals)^2)
l.BETA.st <- rmvt(K, sigma = s2*solve((t(X) %*% X)), df =df, delta = ols.l.beta )
l.S2.st <- rinvgamma(K, df/2, s2/2)


save.image(file="/Users/MacUser/Desktop/MA578/Final Project/analysis.RData")
save(BETA, LAMBDA, S2, BETA.st, l.BETA, l.BETA.st, S2.st, l.S2, l.LAMBDA, l.S2, l.S2.st, file="/Users/MacUser/Desktop/MA578/Final Project/AnalysisMain.RData")


