# PB stochastic simulation of a random walk in continuous time

# Functions dmarkov (cmarkov) generates a sample path of a Markovian random walk 
# in discrete (continuous) time. 
# Arguments:
# S -     an M times N (number of species times number of channels) stoichiometric matrix 
# h -     channel weight/rate function h:Z^M -> R^N
# x0-     initial state from Z^M (the initial time is zero).
# t -     a vector of recording time points at which we want the value of the process to be
#         returned. It must be nonegative and increasing.
# Output: 
# x -     a matrix M times lenght(t), which contains the values of a sample path of
#         at the recording timepoints

# DISCRETE TIME MARKOV CHAIN

dmarkov <- function(S, h, x0, t){
t_cur = 0
x_cur = x0
t_index = 1
M = nrow(S)
N = nrow(S)
x = matrix(nrow=M, ncol=length(t))
while (t_cur <= t[length(t)]) {
    if (t_index <= length(t) & t[t_index] == t_cur) {
		x[,t_index] <- x_cur
		t_index <- t_index + 1
	}
	t_cur <- t_cur + 1
    j = sample(N,1, prob=h(x_cur))
	x_cur = x_cur + S[,j]
}
return(x)
}

dmarkovsample <- function(K, S, h, x0, t){
sam = matrix(nrow=nrow(S), ncol=K)
for (k in 1:K) {
	x <-dmarkov(S, h, x0, c(t))
	sam[,k] <- x[,1]
}
return(sam)
}

# EXAMPLES

# Ehrenfest Discrete

# PART 1 - time trace
t <- seq(100)
x <- dmarkov(matrix(c(-1,1,1,-1), nrow=2, ncol=2), function(x) x, c(20,0), t)
plot(t, x[1,], type="b", main="White balls")
# PART 2 - distribution
K = 10000
sam <- dmarkovsample(K, matrix(c(-1,1,1,-1), nrow=2, ncol=2), function(x) x, c(20,0), 50)
x <- seq(0,20,2) # step by 2 because of periodicity
plot(table(sam[1,]))
lines(x, K*2*dbinom(x, 20, 0.5),type="b") # factor 2 because of periodicity

# Polya-Eggenberger

# PART 1 - time trace
t <- seq(100)
x <- dmarkov(matrix(c(1,0,0,1), nrow=2, ncol=2), function(x) x, c(1,1), t)
plot(t, x[1,], type="b", main="White balls")
# PART 2 - multiple time traces
t <-seq(100)
plot(0,0,xlim=c(0,100),ylim=c(0,101), xlab="Time", ylab="White balls", col="white")
for (i in 1:50) {
x <- dmarkov(matrix(c(1,0,0,1), nrow=2, ncol=2), function(x) x, c(1,1), t)
lines(t,x[1,])
}
# PART 3 - histogram
sam <- dmarkovsample(10000, matrix(c(1,0,0,1), nrow=2, ncol=2), function(x) x, c(1,1), 100)
hist(sam[1,], breaks=seq(1,101,10))

# CONTINUOUS TIME MARKOV CHAIN

cmarkov <- function(S, h, x0, t, iter_max=Inf){
t_cur <- 0.0
x_cur <- x0
iter <- 1
t_index <- 1
M <- nrow(S)
N <- ncol(S)
x <- matrix(nrow=M, ncol=length(t))
while (t_cur < t[length(t)] & iter <= iter_max) {

    rate = sum(h(x_cur))
    if (rate > 0) {
        t_cur <- t_cur + rexp(n=1, rate=rate)
        j <- sample(N,1, prob=h(x_cur))
    } else {
        t_cur <- Inf
        j <- 1
    }

	while (t_index <= length(t) & t[t_index] < t_cur) {
		x[,t_index] <- x_cur
		t_index <- t_index + 1
	}	

	x_cur <- x_cur + S[,j]
	iter <- iter + 1
}
return(x)
}

cmarkovsample <- function(K, S, h, x0, t, iter_max=Inf){
sam = matrix(nrow=nrow(S), ncol=K)
for (k in 1:K) {
	x <-cmarkov(S, h, x0, c(t), iter_max)
	sam[,k] <- x[,1]
}
return(sam)
}

# EXAMPLES

# Ehrenfest Continuous

# PART 1 - time trace
t <- seq(0,10,0.005)
x <- cmarkov(matrix(c(-1,1,1,-1), nrow=2, ncol=2), function(x) x, c(20,0), t)
plot(t, x[1,], type="l", main="White balls")
# PART 2 - distribution
K = 10000
sam <- cmarkovsample(K, matrix(c(-1,1,1,-1), nrow=2, ncol=2), function(x) x, c(20,0), 10)
x <- seq(0,20,1) # step by 1 because of periodicity
plot(table(sam[1,]))
lines(x, K*dbinom(x, 20, 0.5),type="b") # factor 1 because of periodicity

# Lotka Volterra - parameters adopted from Wilkinson's

th <- c(1, 0.005, 0.6)
x0 <- c(50, 100)
S <- matrix(c(1,0,-1,1,0,-1), nrow=2, ncol=3)
h <- function(x) c(th[1]*x[1], th[2]*x[1]*x[2], th[3]*x[2])
t <- seq(0,100,0.1)
x <- cmarkov(S, h, x0, t)
par(mfcol=c(2,2))
plot(t,x[1,], type="l", main="KorisĹĄ")
plot(t,x[2,], type="l", main="Dravec")

library("deSolve")
LVD <- function(t, x, th) list(c(th[1]*x[1] - th[2]*x[1]*x[2], th[2]*x[1]*x[2] - th[3]*x[2]))
res <- ode(x0, t, LVD, th)
#par(mfrow=c(2,1))
plot(res[,1],res[,2], type="l", main="KorisĹĄ")
plot(res[,1],res[,3], type="l", main="Dravec")

# Gene expression

# PART1: time traces
th <- c(1,5,100,1)
x0 <-c(0,0)
S <- matrix(c(1,0,-1,0,0,1,0,-1), nrow=2, ncol=4)
h <- function(x) c(th[1], th[2]*x[1], th[3]*x[1], th[4]*x[2])
t <- seq(0,10,0.005)
x <- cmarkov(S, h, x0, t)
par(mfrow=c(2,1))
plot(t,x[1,], type="l", main="mRNA")
plot(t,x[2,], type="l", main="Protein")
# PART2: moments
sam <- cmarkovsample(1000,S,h,x0,10)
mean(sam[1,])
mean(sam[2,])
sd(sam[1,])
sd(sam[2,])
sd(sam[1,])/mean(sam[1,])
sd(sam[2,])/mean(sam[2,])
sd(sam[1,])**2/mean(sam[1,])
sd(sam[2,])**2/mean(sam[2,])

# SIR

# PART1: definition
R0 = 2.0
#R0 = 0.5
N0 = 1000
th = c(R0/N0, 1)
x0 = c(N0-1, 1)
S <- matrix(c(-1,1,0,-1), nrow=2, ncol=2)
h <- function(x) c(th[1]*x[1]*x[2], th[2]*x[2])
# PART1: time traces
t = seq(0,20,0.05)
x <- cmarkov(S, h, x0, t)
par(mfcol=c(3,2))
plot(t, x[1,], type="l", main="Susceptible")
plot(t, x[2,], type="l", main="Infected")
plot(t, N0 - x[1,] - x[2,], type="l", main="Recovered")

SIRD <- function(t, x, th) list(c(-th[1]*x[1]*x[2], th[1]*x[1]*x[2] - th[2]*x[2]))
res <- ode(x0, t, SIRD, th)
plot(res[,1], res[,2], type="l", main="Susceptible")
plot(res[,1], res[,3], type="l", main="Infected")
plot(res[,1], N0 - res[,2] - res[,3], type="l", main="Recovered")

# PART2: epidemic size
sam <- cmarkovsample(1000, S, h, x0, 20)
R = N0 - sam[1,] - sam[2,]
hist(R, breaks=20)
hist(R[R < 400], breaks=20)
hist(R[R > 400], breaks=20)
