# Stochasticke simulacne metody
# Radoslav Harman, KAMS FMFI UK
# Cvicenie 6 (Monte Carlo vypocet integralov) - strucny zaznam prikazov pre R

# Poznamka: Tento vypocet sluzi na jednoduchu ilustraciu metod; jednorozmerne integraly
# sa obvykle daju pocitat ovela efektivnejsie klasickymi numerickymi metodami.

mc.crude <- function(h, lam, ru, n)
{
  # "Obycajna" Monte-Carlo metoda
  # h ... integrovana funkcia
  # lam ... Lebesguova miera I (mnozina, po ktorej integrujeme h)
  # ru ... generator z rovnomerneho rozdelenia na I
  # n ... pocet realizacii

  Z <- lam * h(ru(n))
  I <- mean(Z)

  list(I=I, delta=qnorm(0.975) * sd(Z) / sqrt(n))
}

h <- function(x) {
  exp(-x^2 / 2)
}


mc.imp <- function(h, f, rf, n) {
  # Monte-Carlo metoda vahoveho vyberu
  # h ... integrovana funkcia
  # f ... predpis simulacnej hustoty
  # rf ... generator zo simulacnej hustoty
  # n ... pocet realizacii

  X <- rf(n)
  Z <- h(X)/f(X)
  I <- mean(Z)

  list(I=I, delta=qnorm(0.975) * sd(Z) / sqrt(n))
}

dexpc <- function(x) {
  # Hustota rozdelenia Exp(1/2) useknuteho na (0,1)
  (1/2) * exp(-x/2) / (1 - exp(-1/2))
}

rexpc <- function(n) {
  # Generator z rozdelenia Exp(1/2) useknuteho na (0,1)
  -2 * log(1- (1 - exp(-0.5)) * runif(n))
}

system.time(print(mc.crude(h, 1, runif, 1000000)))
system.time(print(mc.imp(h, dexpc, rexpc, 1000000)))

integrate(h, 0, 1)
(pnorm(1) - 0.5) * sqrt(2 * pi)


# Monte Carlo odhad cisla pi

# Najjednoduchsia MC metoda konzistentneho odhadu cisla pi 

est.pi.crude <- function (n) 
{
  # Najjednoduchsia metoda MC odhadu cisla pi
  4 * sum(runif(n)^2+runif(n)^2<1) / n
}

est.pi.crude(10)
est.pi.crude(1000)
est.pi.crude(100000)
est.pi.crude(10000000)
pi
# Du: najdite metodu vypoctu pribliznych intervalov spolahlivosti

# Specialna variance reduction metoda z knihy Ross: Simulation 
est.pi <- function(n) 
{
  # Monte-Carlo estimation of pi
  U <- runif(n)
  2 * mean(sqrt(1-((U+(1:n)-1)/n)^2)+sqrt(1-(((1:n)-U)/n)^2))
}
 
options(digits=20)
est.pi(10)
est.pi(1000)
est.pi(100000)
est.pi(10000000)
pi
# Du: najdite metodu vypoctu pribliznych intervalov spolahlivosti