# Stochasticke simulacne metody
# Radoslav Harman, KAMS FMFI UK
# Cvicenie 3 (Generovanie spojitych premennych a vektorov) - strucny zaznam prikazov pre R

# Poznamka: Nasledovne funkcie sluzia na demonstraciu teoretickych principov generovania.
# Na prakticke simulacie v jazyku R je vyhodnejsie pouzit standardne funkcie rnorm, rbeta a rgamma.


rnorm.rej <- function(n, mu=0, sigma=1) {
   # Specialnou zamietacou metodou generuje n iid realizacii z N(mu,sigma^2)
   
   res <- rep(0, n)
   for (i in 1:n) {
      ok <- FALSE
      while (!ok) {
         Y <- rexp(1)
         if (runif(1) < exp(-(Y - 1)^2/2)) ok <- TRUE
      }
      res[i] <- Y
   }
   rsgn <- 2*floor(2*runif(n)) - 1
   
   return(sigma * rsgn * res + mu)
}

res <- rnorm.rej(10000); hist(res); ks.test(res, pnorm)


rnorm.bm <- function(n, mu=0, sigma=1) {
   # Boxovou-Muellerovou metodou generuje n iid realizacii z N(mu,sigma^2)
   
   m <- floor(n/2) + 1; res <- c(0, 2*m)
   for (i in 1:m) {
      U <- runif(1); L <- rexp(1)
      res[2*i - 1] <- sqrt(2*L)*cos(2*pi*U)
      res[2*i] <- sqrt(2*L)*sin(2*pi*U)
   }
   
   return(sigma*res + mu)
}

res <- rnorm.bm(10000); hist(res); ks.test(res, pnorm)


rbeta.rej <- function(n, a=1, b=1) {
   # Von-Johnkovou metodou generuje n iid realizacii z Be(a,b)
   
   res <- rep(0, n)
   for (i in 1:n) {
      ok <- FALSE
      while (!ok) {
         X <- runif(2)
         V <- X[1]^(1/a); W <- X[2]^(1/b)
         if (V + W < 1) ok <- TRUE
      }
      res[i] <- V/(V + W)
   }
   
   return(res)
}

res <- rbeta.rej(10000, a = 2, b = 2); hist(res)
ks.test(res, pbeta, shape1 = 2, shape2 = 2)


rgamma.bex <- function(n, a=1, b=1)
{
   # Generuje n iid realizacii z Gamma(a,b)
   
   m <- floor(b); r <- b - m
   res <- rep(0, n)
   for (i in 1:n) {
      Y <- 0; if (r > 0) Y <- rbeta(1, shape1 = r, shape2 = 1 - r)
      Z <- 0; if (m > 0) Z <- sum(rexp(m))
      res[i] <- (Z + rexp(1)*Y)/a
   }
   
   return(res)
}

res <- rgamma.bex(10000, a = 2, b = 3); hist(res)
ks.test(res, pgamma, rate = 2, shape = 3)


rnorm.2d <- function(n, mu1=0, mu2=0, sd1=1, sd2=1, rho=0) {
   # Generátor realizacii vektora (X1,X2) zo zdruzene normalneho rozdelenia
   # mu1=E(X1), mu2=E(X2), sd1^2=D(X1), sd2^2=D(X2), rho=cor(X1,X2)
   
   X1 <- rnorm(n); X2 <- rnorm(n)
   Y1 <- sd1*X1 + mu1
   Y2 <- sd2*(rho*X1 + sqrt(1 - rho^2)*X2) + mu2
   
   return(rbind(Y1, Y2))
}

res <- rnorm.2d(5000); plot(res[1,], res[2,], pch = 19, cex = 0.1)
res <- rnorm.2d(5000, rho = 0.5); plot(res[1,], res[2,], pch = 19, cex = 0.1)
res <- rnorm.2d(5000, rho = -0.9); plot(res[1,], res[2,], pch = 19, cex = 0.1)


rsimplex <- function(n, m) {
   # Generátor n realizacii z rovnomerneho rozdelenia na m-rozmernom simplexe
   # Poznamka: klasicky, nie pravdepodobnostny simplex
   
   res <- matrix(ncol = n, nrow = m)
   for (i in 1:n) {
      X <- runif(m); Y <- sort(X)
      Z <- c(Y[1], Y[2:m] - Y[1:(m - 1)])
      res[, i] <- Z
   }
   
   return(res)
}

res <- rsimplex(10000, 2); plot(res[1,], res[2,], pch = 19, cex = 0.1, asp = 1) # +m=3,10


rsphere <- function(n, m) {
   # Generátor n realizacii z rovnomerného rozdelenia na povrchu m-rozmernej gule
   
   res <- matrix(ncol = n, nrow = m)
   for (i in 1:n) {
      X <- rnorm(m)
      res[, i] <- X/sqrt(sum(X^2))
   }
   
   return(res)
}

res <- rsphere(10000, 3); plot(res[1,], res[2,], pch = 19, cex = 0.1, asp = 1) # +m=4,10


rball <- function(n, m) {
   # Generátor n realizacii z rovnomerného rozdelenia vo vnutri m-rozmernej gule
   
   res <- matrix(ncol = n, nrow = m)
   for (i in 1:n) {
      X <- rnorm(m)
      res[,i] <- runif(1)^(1/m)*X/sqrt(sum(X^2))
   }
   
   return(res)
}

res <- rball(10000, 2); plot(res[1,], res[2,], pch = 19, cex = 0.1, asp = 1) # +m=3,4,10


rtorus <- function(n, R, r) {
   # Evkin generator n rovnomernych bodov na hranici torusu s polomermi R, r
   
   rphitor <- function(n, R, r) {
      # Generuje n hodnot s rozdelenim uhla phi pre torus s polomermi R, r
      res <- rep(0, n); k <- 0; OK <- FALSE
      while (!OK) {
         U <- 2*pi*runif(1); V <- (R + r)*runif(1)
         if (abs(R + r*cos(U)) > V) {
            k <- k + 1; if (k == n) OK <- TRUE
            res[k] <- U
         }
      }
      return(res)
   }
   
   phi <- rphitor(n, R, r)
   psi <- 2*pi*runif(n)
   x <- cos(psi)*(R + r*cos(phi))
   y <- sin(psi)*(R + r*cos(phi))
   z <- r*sin(phi)
   
   return(rbind(x, y, z))
}

H <- matrix(rnorm(9), ncol = 3); U <- eigen(H %*% t(H))$vectors
plot(t((U %*% rtorus(10000, 5, 1))[1:2,]), pch = 19, cex = 0.2, asp = 1)
