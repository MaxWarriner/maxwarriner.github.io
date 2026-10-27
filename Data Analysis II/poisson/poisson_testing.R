dat <- ISLR2::Bikeshare
X <- as.matrix(dat[,9:10])
y <- dat[, 15]

tolerance <- 1e-6      # Convergence threshold
#y <- # response
offset <- log(1) # 1 is no offset
#X <- # design matrix
betas <- matrix(0, nrow = ncol(X), ncol = 1) # Initialize coefficients to 0
ns <- 1 # Initialize n
nbinom.ll <- function(ns, y, mu){-sum(dnbinom(y, size=ns, mu=mu, log=T))}
repeat{
  eta <- offset + X %*% betas
  mu.hat.i <- exp(eta)
  ns <- optim(par=ns, fn=nbinom.ll, method="Brent", lower=0, upper=10000, y=y, mu=mu.hat.i)$par
  W <- diag(as.vector(mu.hat.i^2 / (1/ns * mu.hat.i^2 + mu.hat.i)))
  z <- eta + (y - mu.hat.i) / (1/ns * mu.hat.i^2 + mu.hat.i)
  Hessian.inv <- solve(t(X) %*% W %*% X)
  betas.new <- Hessian.inv %*% t(X) %*% W %*% z
  # Check for Convergence
  diff <- sum((betas.new - betas)^2) # sums of squared differences
  betas <- betas.new
  std.errors <-sqrt(diag(Hessian.inv))
  if (diff < tolerance) {
    break
  }
}
cbind(betas, std.errors)
summary(glm(bikers ~ temp + atemp, family = poisson, data = dat))



data(Affairs, package = 'AER')
X <- as.matrix(Affairs[,3:4])
y <- Affairs[, 1]

tolerance <- 1e-6      # Convergence threshold
#y <- # response
offset <- log(1) # 1 is no offset
#X <- # design matrix
betas <- matrix(0, nrow = ncol(X), ncol = 1) # Initialize coefficients to 0
ns <- 1 # Initialize n
nbinom.ll <- function(ns, y, mu){-sum(dnbinom(y, size=ns, mu=mu, log=T))}
repeat{
  eta <- offset + X %*% betas
  mu.hat.i <- exp(eta)
  ns <- optim(par=ns, fn=nbinom.ll, method="Brent", lower=0, upper=10000, y=y, mu=mu.hat.i)$par
  W <- diag(as.vector(mu.hat.i^2 / (1/ns * mu.hat.i^2 + mu.hat.i)))
  z <- eta + (y - mu.hat.i) / (1/ns * mu.hat.i^2 + mu.hat.i)
  Hessian.inv <- solve(t(X) %*% W %*% X)
  betas.new <- Hessian.inv %*% t(X) %*% W %*% z
  # Check for Convergence
  diff <- sum((betas.new - betas)^2) # sums of squared differences
  betas <- betas.new
  std.errors <-sqrt(diag(Hessian.inv))
  if (diff < tolerance) {
    break
  }
}
cbind(betas, std.errors)
summary(glm(affairs ~ age + yearsmarried, family = poisson, data = Affairs))





