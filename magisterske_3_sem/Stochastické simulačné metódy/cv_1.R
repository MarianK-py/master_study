# Stochasticke simulacne metody
# Radoslav Harman, KAMS FMFI UK
# Cvicenie 1 (Testy RNG) - strucny zaznam prikazov pre R

# Du: pozrite si ako vyzeraju bezne dostupne hardwarove generatory nahodnych cisel, napr:
# https://www.amazon.com/TrueRNGpro-Hardware-Random-Number-Generator/dp/B01JTJ6D0S/
# https://www.amazon.com/TrueRNG-V3-Hardware-Random-Generator/dp/B01KR2JHTA/
# Vyhladajte si akou rychlostou generuju.

# Du : Pozrite si tiez stranky www.random.org a https://qrng.anu.edu.au/

# PB: nic z toho som nekukal

RNGkind()                   # PB: toto zisti aky generator fici v R
.Random.seed                # PB: toto funguje az po nastaveni seedu
set.seed(123456789)         # Velmi dolezite pre replikaciu a debugging. Nepouzivat v cykle!

# PB: Reset seed
set.seed(NULL)
set.seed(Sys.time())

# Du: Precitajte si podrobnejsie https://stat.ethz.ch/R-manual/R-patched/library/base/html/Random.html

# PB: du som ignoroval

lcg <- function(a, c, m, x0, n) {
   # Elementarna implementacia LCG prveho radu
   # a, c, m, x0 ... multiplikator, inkrement, modulus, seed
   # n ... dlzka postupnosti pseudonahodnych cisiel
   
   x <- c(x0, rep(0, n - 1))
   for (i in 2:n) x[i] <- (a * x[i - 1] + c) %% m
   
   return(x)
}

x0 <- lcg(13, 0, 32, 17, 33); print(x0)                    # Len na "overenie" teoretickej vety z prednasky (PB: Veta 2.3, prva cast - nulovy intercept)
x0 <- (17*13**(seq(1,32)-1)) %% 32; print(x0)              # PB: teoreticky staci modulo zobrat nakoniec; prakticky to nefunguje
x0 <- lcg(17, 17, 32, 17, 33); print(x0)                   # PB: Veta 2.3, druha cast, nenulovy intercept (za Tabulkou 1)


# PB: Tu zacina Sekcia 2.5

x1 <- lcg(17, 123456, 2^31 - 1, 123456, 500*2^12) / (2^31 - 1) # LCG Generator s "nahodne" zvolenymi parametrami
x2 <- lcg(7^5, 0, 2^31 - 1, 123456, 500*2^12) / (2^31 - 1)     # Generator "Minimal standard"
x3 <- runif(500 * 2^12)                                        # Generator Rka

# Vizualna inspekcia histogramu
par(mfrow = c(3, 1))
hist(x1, ylim = c(0.95, 1.05), freq = FALSE); lines(c(0, 1), c(1, 1), col = "red", lty = "dashed")
hist(x2, ylim = c(0.95, 1.05), freq = FALSE); lines(c(0, 1), c(1, 1), col = "red", lty = "dashed")
hist(x3, ylim = c(0.95, 1.05), freq = FALSE); lines(c(0, 1), c(1, 1), col = "red", lty = "dashed")

# Vizualna inspekcia realizacii dvojrozmernych subvektrov 
# PB: toto je device-dependent. pod linux funguje x11()
windows(); plot(x1[1:(250 * 2^8)], x1[(250 * 2^8 + 1):(500 * 2^8)], pch = 19, cex = 0.1)
windows(); plot(x2[1:(250 * 2^8)], x2[(250 * 2^8 + 1):(500 * 2^8)], pch = 19, cex = 0.1)
windows(); plot(x3[1:(250 * 2^8)], x3[(250 * 2^8 + 1):(500 * 2^8)], pch = 19, cex = 0.1)

# Kolmogorov-Smirnov test
ks.test(x1, punif)
ks.test(x2, punif)
ks.test(x3, punif)

# Jednorozmerny test dobrej zhody
chisq.test(table(ceiling(1000 * x1)), p = rep(1/1000, 1000))
chisq.test(table(ceiling(1000 * x2)), p = rep(1/1000, 1000))
chisq.test(table(ceiling(1000 * x3)), p = rep(1/1000, 1000))

# PB - vynechat poker test? skripta hovoria, ze moc pouzivany neni

poker.test <- function(x) {
   # Takzvany trojrozmerny pokrovy test
   # (Specialny test dobrej zhody s rovnomernym rozdelenim na [0,1]^m, m=3.)
   # x ... vektor testovanych cisel v intervale [0,1]
   
   d <- length(x); d <- d - (d %% 3)
   x <- x[1:d]; n <- d/3
   
   N <- c(0, 0, 0)
   for (l in 1:n) {
      s <- length(unique(floor(3 * x[(3 * l - 2):(3 * l)])))
      N[s] <- N[s] + 1
   }
   
   # V tomto teste mame k=3 oblasti a nasledovne p su p_1,p_2,p_3 z prednasky
   return(chisq.test(N, p = c(3/27, 18/27, 6/27)))
}

# Trojrozmerny pokrovy test
# PB - x1 fails but only marginally - vynechat
poker.test(x1) 
poker.test(x2)
poker.test(x3)

# PB - serialny test robit. x1 fails grandly.

serial.test <- function(x, k) {
   # Test nahodnosti pomocou serialnej korelacie
   # x ... testovany vektor
   # k ... "posun"
   
   n <- length(x)
   rk <- cor(x[(k + 1):n], x[1:(n - k)])
   print(paste("Koeficient serialnej korelacie:", rk), quote = FALSE)
   pval <- 2 * (1 - pnorm(abs(sqrt(n) * rk)))
   
   print(paste("p-hodnota testu:", pval), quote = FALSE)
}

# Test serialnej korelacie
serial.test(x1, 2) 
serial.test(x2, 2)
serial.test(x3, 2)

# PB: Sposob ako sa presvedcit o "mriezkovej strukture" (Priklad 2.4) x1 fails

> X1 = matrix(x1[1:(500*2**8)], ncol=2, byrow=TRUE)
> x11(); plot(X1[,1],X1[,2], pch=19, cex=0.1)

> X2 = matrix(x2[1:(500*2**8)], ncol=2, byrow=TRUE)
> x11(); plot(X2[,1],X2[,2], pch=19, cex=0.1)

# Idea ilustrovana na hracej matici
> X1 = matrix(seq(1,10), ncol=2, byrow=TRUE)
> X1
     [,1] [,2]
[1,]    1    2
[2,]    3    4
[3,]    5    6
[4,]    7    8
[5,]    9   10
> X1[,1]
[1] 1 3 5 7 9
> X1[,2]
[1]  2  4  6  8 10

# PB - robit dalsie testy? Perhaps birthday test because the 'minimal standard' fails it.

extreme.test <- function(x) {
   # Test nezavislosti pomocou extremalnych bodov
   #  x ... testovany vektor cisel
   
   n <- length(x)
   y <- 0
   for (i in 2:(n - 1)) {
      ismax <- (x[i - 1] < x[i]) & (x[i] > x[i + 1])
      ismin <- (x[i - 1] > x[i]) & (x[i] < x[i + 1])
      if (ismax | ismin) y <- y + 1
   }
   print(paste("Pocet extremalnych bodov v postupnosti:", y), quote = FALSE)
   
   # ak je x generovane z idealneho RNG, tak y ma asymptoticky normalne rozdelenie
   # s nasledovnou strednou hodnotou a disperziou
   Ey <- 2 * (n - 2) / 3; Dy <- (16 * n - 29) / 90
   
   pval <- 2 * (1 - pnorm(abs((y - Ey) / sqrt(Dy))))
   print(paste("p-hodnota testu:", pval), quote = FALSE)
}

# Test extremalnych hodnot 
# PB - x1 fails grandly but that already failed previous tests - vynechat
extreme.test(x1) 
extreme.test(x2)
extreme.test(x3)


birthdays.test <- function(x) {
   # Specifikacia testu "birthday spacings"
   # x ... testovany vektor cisel z [0,1]
   # x v tejto implementacii musi mat dlzku 500*2^12
   # Pocet "rokov" je N=500, "dni" v roku je d=2^32, "narodenin/ludi" v roku je k=2^12
   # Za kazdy rok sa vypocita vektor odstupov (v dnoch) medzi po sebe nasledujucimi narodeninami
   # Potom sa urci hodnota kolko tychto odstupov sa zopakovalo (viac ako raz)
   # Pre korektny generator je rozdelenie tychto hodnot priblizne Po(k^3/(4*d))=Po(4)
   
   N <- 500; k <- 2^12; d <- 2^32
   res <- rep(0, N)
   for (i in 1:N) {
      z <- sort(ceiling(d * x[((i - 1) * k + 1):(i * k)]))
      sp <- table(z[2:k] - z[1:(k - 1)])
      res[i] <- length(sp[sp > 1])
   }
   
   nn <- rep(0, 8)
   for (j in 0:7) nn[j + 1] <- length(res[res == j])
   nn <- c(nn, length(res[res > 7]))
   pp <- dpois(0:7, lambda = 4); pp <- c(pp, 1 - sum(pp))
   
   return(chisq.test(nn, p = pp))
}

# Birthday spacings test
# PB: x1 and x2 fail grandly
birthdays.test(x1)
birthdays.test(x2)
birthdays.test(x3)

# Ovela dokladnejsie testovanie je mozno vykonat pomocou balikov:
# https://en.wikipedia.org/wiki/Diehard_tests
# https://en.wikipedia.org/wiki/TestU01