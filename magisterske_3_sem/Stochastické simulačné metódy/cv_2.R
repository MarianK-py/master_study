# Stochasticke simulacne metody
# Radoslav Harman, KAMS FMFI UK
# Cvicenie 2 (Generovanie diskretnych premennych a vektorov) - strucny zaznam prikazov pre R

rdisc.cut <- function(n, p, d) {
   # Metoda deliacich bodov (angl. cutpoints method)
   # Generuje n realizacii nahodnej premennej nadobudajucej
   # hodnoty 1,...,n s pravdepodobnostami p[1],...,p[m]
   # d ... pocet deliacich bodov (parameter metody)
   
   # Vypocet kumulativnych pravdepodobnosti
   m <- length(p)
   s <- c(p[1], rep(0, m - 1))
   for (i in 2:m) s[i] <- s[i - 1] + p[i]
   
   # Vypocet deliacich bodov
   I <- rep(0, d)
   for (j in 1:d) I[j] <- min((1:m)[s > (j - 1) / d])
   
   # Generovanie pomocou predvypocitanych deliacich bodov
   res <- rep(0, n)
   for (i in 1:n) {
      U <- runif(1)
      k <- I[ceiling(d * U)]
      while (s[k] <= U) k <- k + 1
      res[i] <- k
   }
   
   return(res)
}

res <- rdisc.cut(1000, rep(1/10, 10), 3) 
hist(res, breaks = 100)
chisq.test(table(res), p = rep(1/10, 10))

pb <- dbinom(0:1000, size = 1000, prob = 0.5); plot(pb, type = "h")
system.time(invisible(rdisc.cut(100000, pb, 1)))
system.time(invisible(rdisc.cut(100000, pb, 10)))
system.time(invisible(sample(0:1000, 100000, replace = TRUE, prob = pb)))
system.time(invisible(rbinom(100000, size = 1000, prob = 0.5)))

rwor <- function(n, m) {
   # Generuje rovnomerny nahodny vyber n prvkov bez navratu z mnoziny {1,...,m}
   # (Da sa chapat ako generator z rozdelenia specialneho diskretneho nahodneho vektora.)
   # Pre n=m generuje rovnomernu nahodnu permutaciu cisel 1,...,m.
   
   r <- 1:m; res <- rep(0, n)
   for (i in 1:n) {
      l <- m - i + 1; k <- sample(1:l, 1)
      res[i] <- r[k]; if (k < l) r[k] <- r[l]
   }
   
   return(res)
}

rwor(6, 49); sample(1:49, 6)
rwor(49, 49); sample(1:49, 49)

rwor.reservoir <- function(n, m) {
   # Generuje rovnomerny nahodny vyber n prvkov bez navratu z mnoziny {1,...,m}
   # Pouziva algoritmus s nazvom reservoir sampling
   
   res <- 1:n
   for (i in (n + 1):m) {
      if (runif(1) < n/i)  
         res[sample(1:n, 1)] <- i
   }
   
   return(res[sample(1:n)])
}

rwor.reservoir(6, 49)

rmult <- function(n, N, p) {
   # Generuje n realizacii z multinomickeho rozdelenia Mult(N,p[1],...,p[m])
   
   m <- length(p)
   res <- matrix(ncol = n, nrow = m)
   for (i in 1:n) {
      X <- rep(0, m); X[1] <- rbinom(1, size = N, prob = p[1])
      Nrem <- N - X[1]; prem <- 1 - p[1]
      for (j in 2:m) {
         X[j] <- rbinom(1, size = Nrem, prob = min(p[j] / prem, 1))
         Nrem <- Nrem - X[j]; prem <- prem - p[j]
      }
      res[, i] <- X
   }
   
   return(res)
}

rmult(10, 10, c(0.1, 0.2, 0.3, 0.4))
