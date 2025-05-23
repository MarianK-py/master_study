library(astsa)
library(urca)
library(areaplot)
library(fGarch)
library(polynom)

# Kratke funkcie na lahsie vykreslovnanie
fn_plot_train_and_pred <- function(train, predict, model) {
  plot(as.integer(rownames(train)), train$interestRate,  col="black", xlim = c(90, 370), ylim = c(- 4, 6))
  lines(as.integer(rownames(train)), data_train$interestRate, col="black")
  
  confplot(as.integer(rownames(predict)), model$pred-2*model$se, model$pred+2*model$se, add=TRUE)
  confplot(as.integer(rownames(predict)), model$pred-model$se, model$pred+model$se, col="darkgray", add=TRUE)
  
  points(as.integer(rownames(predict)), predict$interestRate,  col="green")
  lines(as.integer(rownames(predict)), predict$interestRate, col="green")
  points(as.integer(rownames(predict)), model$pred, col="red")
  lines(as.integer(rownames(predict)), model$pred, col="red")
  
}
fn_plot_pred <- function(predict, model) {
  confplot(as.integer(rownames(predict)), model$pred-2*model$se, model$pred+2*model$se, ylim = c(- 4, 4))
  confplot(as.integer(rownames(predict)), model$pred-model$se, model$pred+model$se, col="darkgray", add=TRUE)
  
  points(as.integer(rownames(predict)), predict$interestRate,  col="green")
  lines(as.integer(rownames(predict)), predict$interestRate, col="green")
  points(as.integer(rownames(predict)), model$pred, col="red")
  lines(as.integer(rownames(predict)), model$pred, col="red")
  
}


# Ako prve si nastavime cestu k nasim datam a nacitame ich
setwd("C:/Users/majko/Desktop/Casove_rady_projekt")
data <- read.csv("ECB Data Portal_20231112202535.csv",
                 header=TRUE) 
# Teraz si tie data pripravime (pridame stlpec roku)
data <- data[c(2,3)]
colnames(data) <- c("yearMonth", "interestRate")
data["year"] <- as.integer(substring(data$yearMonth, 1,4 ))
# Tieto data si to
data_train <- data[data$year %in% 2001:2020,]
data_predict <- data[data$year > 2020,]

# ako prve si vykreslime nase trenovacie data
plot(as.integer(rownames(data_train)), data_train$interestRate)
lines(as.integer(rownames(data_train)),data_train$interestRate)
# vidime klesajuci trend v datach ktory si este skusme potvrdit ACF funciou
acf2(data_train$interestRate)

# ACF klesa velmi pomaly co potvrdzuje trend preto bude pracovat dalej uz len 
# s diferenciami nasich dat, tie si takisto vykrslime
l = length(data_train[,1])
plot(as.integer(rownames(data_train))[1:(l-1)],diff(data_train$interestRate))
lines(as.integer(rownames(data_train))[1:(l-1)],diff(data_train$interestRate))
acf2(diff(data_train$interestRate))
# vidime, ze ani data ani ACF nevykazuje znamky trendu preto nemusi 
# pocitat dalsiu diferenciu

# Dalsi krok je testovani jednotkovy koren pomocou ADF testu
# v tomto pripade pouzijeme type="drift" kedze tvrdime, ze nase povodne
# data maju linearny trend a preto nase diferencie budu mat konstantny clen
urTest = ur.df(diff(data_train$interestRate), type = "drift", lags = 12, selectlags = "BIC")
summary(urTest)
# Vidime, ze hodnota statisiky je vyrazne mensia ako 5% kriticka hodnota
# preto mozeme nulovu hypotezu, ze existuje jednotkovy koren zamietnut
# a tym padom nemusi nase data dalej diferencovat

Box.test(diff(data_train$interestRate),
         lag = 12,
         type = "Ljung-Box") 

testik <- rnorm(100)^5

lb <- mapply(function(x) Box.test(diff(data_train$interestRate), lag = x, type = "Ljung-Box")$p.value , 1:12)
lb <- mapply(function(x) Box.test(testik , lag = x, type = "Ljung-Box")$p.value , 1:12)


plot(lb, ylim=c(0,1))
lines(1:12,rep(0.05, 12), lty='dashed', col = "blue")

# Teraz sa mozeme pustit do hladania najlepsieho modelu
acf2(diff(data_train$interestRate))

sarima(data_train$interestRate, 1,1,0, no.constant = FALSE)
sarima(data_train$interestRate, 0,1,1, no.constant = FALSE)


sarima(data_train$interestRate, 1,1,1, no.constant = FALSE)


pred = sarima.for(data_train$interestRate, n.ahead = 33, 0, 1, 1)
lines(as.integer(rownames(data_predict)), data_predict$interestRate, col="green")

fn_plot_train_and_pred(data_train, data_predict, pred)

fn_plot_pred(data_predict, pred)


sarima(data_train$interestRate, 2,0,0, no.constant = FALSE)
pred = sarima.for(data_train$interestRate, n.ahead = 33, 2, 0, 0)




# skuska skusky

library(astsa)
y <- ts(fmri$L1T2[, 2])
plot(y)

acf2(y)

urTest = ur.df(y, type = "drift", lags = 30, selectlags = "BIC")
summary(urTest)

sarima(y, 2,0,0, no.constant = FALSE)



library(astsa)
x <- ts(fmri$L2T2[, 2])

plot(x)

lb <- mapply(function(z) Box.test(x, lag = z, type = "Ljung-Box")$p.value , 1:24)

plot(lb, ylim=c(0,1), xlab = "Lag", ylab = "P-hodnota")
lines(1:24,rep(0.05, 24), lty='dashed', col = "blue")

acf(x)[1]

sarima(x, 0,0,3, no.constant = FALSE)

sarima.for(x, n.ahead = 1, 0, 0, 3)$pred



x = rep(0, 3000)
r = rnorm(3000)

for (i in 2:3000) {
  x[i] = 0.7113787*x[i-1] + r[i]
}

yt <- arima.sim(list(order=c(2,0,0), ar=c(0.5,-0.1), xmean=10), n=5000000)
pacf(yt)[1:2]

yt <- arima.sim(list(order=c(4,0,0), ar=c(0,0,0,0.5)), n=5000000)
print(acf(yt))

ARMAtoMA(ar = c(0.5, 0.2), lag.max=5)


solve(polynomial(coef = c(-1, 1,0,1)))


ww <- garchFit(~garch(1,1), x)


yt <- arima.sim(list(order=c(2,0,3), ar=c(0.01,0.01), ma=c(0.01,0.01,0.01)), n=5000000)
acf(yt)[1:10]
pacf(yt)[1:10]


ARMAacf(ar=c(0.01,0.01), ma=c(0.01,0.01,0.01), lag.max = 5, pacf = TRUE)


library(datasets)
y <- ts(faithful[,2])
plot(y)

mean(y)

urTest = ur.df(y, type = "drift", lags = 30, selectlags = "BIC")
summary(urTest)


acf2(y)

sarima(y, 1,0,0, no.constant = FALSE)


x <- log(AirPassengers)

plot(x)

HoltWinters(x, seasonal = "additive")

h <- HoltWinters(x, alpha=1/3, beta=1/3, gamma=1/3, seasonal = "additive")


vv <- acf(diff(x))


lb <- mapply(function(z) Box.test(diff(x), lag = z, type = "Ljung-Box")$p.value , 1:24)

plot(lb, ylim=c(0,1), xlab = "Lag", ylab = "P-hodnota")
lines(1:24,rep(0.05, 24), lty='dashed', col = "blue")
2024-202.4



yt <- arima.sim(list(order=c(0,0,6), ma=c(-0.2,-0.3,0,-0.24,0.2*0.24,0.3*0.24)), n=5000000)
acf(yt)[1:10]
pacf(yt)[1:10]


yt <- arima.sim(list(order=c(0,0,6), ma=c(-0.2,-0.3,0,-0.24,0.2*0.24,0.3*0.24)), n=5000000)
acf(yt)[1:10]
pacf(yt)[1:10]

yt <- arima.sim(list(order=c(0,0,12), ma=c(rep(0,11),0.1)), n=5000000)
acf(yt)[1:13]
pacf(yt)[1:10]

yt <- arima.sim(list(order=c(0,0,1), ma=c(-0.5549581)), n=5000000)
acf(yt)[1:13]
pacf(yt)[1:10]


solve(polynomial(coef = c(1, 0,0,0, 0.0000001, 0, 0, 0, 0.0000001)))
