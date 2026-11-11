library(pscl)
dat <- pscl::bioChemists
X <- as.matrix(dat[, 4:6])
y <- dat$art

tolerance <- 1e-6      # Convergence threshold
offset <- log(1) # 1 is no offset
X.logit <- X# design matrix
X.pois <- X# design matrix
lambda.hat.i <- exp(offset) # Initialize lambdas
alphas <- matrix(0, nrow = ncol(X.logit), ncol = 1) # Initialize coefficients to 0
betas <- matrix(0, nrow = ncol(X.pois), ncol = 1) # Initialize coefficients to 0
repeat{
  # Logistic Part
  eta.logit <- offset + X.logit %*% alphas 
  pi.hat.i <- exp(eta.logit)/(1+exp(eta.logit))
  delta.i <- (y == 0) * pi.hat.i/(pi.hat.i + (1-pi.hat.i)*exp(-lambda.hat.i)) 
  W.logit <- diag(as.vector(pi.hat.i*(1-pi.hat.i)))
  z.logit <- eta.logit + (delta.i - pi.hat.i) / (pi.hat.i*(1-pi.hat.i))
  Hessian.inv.alphas <- solve(t(X.logit) %*% W.logit %*% X.logit)
  alphas.new <- Hessian.inv.alphas %*% t(X.logit) %*% W.logit %*% z.logit
  # Poisson Part
  eta.pois <- offset + X.pois %*% betas
  lambda.hat.i <- exp(eta.pois)
  W.pois <- diag(as.vector((1-delta.i)*lambda.hat.i))
  z.pois <- eta.pois + (y - lambda.hat.i) / lambda.hat.i
  Hessian.inv.betas <- solve(t(X.pois) %*% W.pois %*% X.pois)
  betas.new <- Hessian.inv.betas %*% t(X.pois) %*% W.pois %*% z.pois
  # Check for Convergence
  diff <- max(sum((alphas.new - alphas)^2), sum((betas.new - betas)^2)) 
  alphas <- alphas.new
  alphas.std.errors <-sqrt(diag(Hessian.inv.alphas))
  betas <- betas.new
  betas.std.errors <-sqrt(diag(Hessian.inv.betas))
  if (diff < tolerance) {
    break
  }
}
cbind(alphas, alphas.std.errors)
cbind(betas, betas.std.errors)

summary(zeroinfl(art ~ kid5 + phd + ment | kid5 + phd + ment, data = dat))


