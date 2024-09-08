setwd("~/Desktop/ThesisData")
data1=read.csv('NewFinalProject.csv')
data1 = read.csv("/home/eoahenkan/Desktop/ThesisData/NewFinalProject.csv", header = TRUE)
View(data1)

## Packages to use for my Research 
library(timeSeries)
library(timeDate)
library(fGarch)
library(zoo)
library(imputeTS)
library(xts)
library(tseries)
library(rugarch)
#library(AER)

library(highcharter)
library(ggplot2)



### PLOTTING RAW DATA


#plot(data)

########################################################################################
## PLOTTING RAW DATA


########################################################################################


#plot(data1 , col= "red" , main = "Plot of the Stock Prices of Bank of Kigali")
#plot(data1[,3]~data1[,1], data=data1, col= "red",  title = "Plot of Bank of Kigali Stock Return", type = "l")

T <- timeSeries(data1[,3], charvec=as.character(data1[,1]), units = NULL, format = NULL, title = "Plot of Bank of Kigali Stock Return")
T
summary(T)
View(T)  
plot(T, col= "blue" , main = " Bank of Kigali Adjusted Closing Prices", xlab="Date",ylab= "Adjusted Closing Price" )


###############################################################################################
## We now apply the imputation method called the Kalman Filter on our Data

###############################################################################################

NewT <- na.kalman(data1[,3], model = "StructTS", smooth = TRUE, nit = -1)
View(NewT)
NewData <- data.frame(data1[,1], NewT)
View(NewData)
#cars <- NewData[1:100,]
#cars
#trydata <- t(cars)

T <- timeSeries(NewData[,2], charvec=as.character(NewData[,1]), units = NULL, format = NULL, title = "Plot of Bank of Kigali Stock Return")
plot(T, col= "red", main = "Plot of the Stock of Adjusted closing Prices", xlab="Date", ylab= "Adjusted Closing Price")
hist(T, col= "red", main = "Histogram of  Adjusted Closing Prices", xlab="Date", ylab= "Adjusted Closing Price")

## Distribution of Miising stock

plotNA.distribution(data1[,3], colPoints = "steelblue", colBackgroundMV = "indianred2", main = "Distribution of Missing Stock Prices of Bank of Kigali", xlab = "Index", ylab = "Adjusted Closing Price", pch = 20, cexPoints = 0.8,
                    col = "black")

## Distribution Bar of Missing Stock Prices

plotNA.distributionBar(data1[,3], breaks = nclass.Sturges(data[,3]), breaksize = NULL, percentage = TRUE, legend = TRUE, axis = TRUE, space = 0,col = c("indianred2", "green2"), main = "Distribution Bar Chart of Missing Stock Prices",
                       xlab = "Time Lapse", ylab = "Percentage")



## New Data
#NewData <- data.frame(data1[,1], NewT)
#NewData

###################################################################################################
## Expected shortfall derived from GARCH(1, 1) models

###################################################################################################

SRLOSS <- timeSeries(-1.0*diff(log(NewData))*100, char.vec = time(NewData))
SRLOSS
## Function for ES of t −GARCH

ESgarch <- function(y, p=0.99){ 
  gfit <- garchFit(formula=~garch(1,1), data = y, cond.dist = "std", trace = FALSE)
  sigma <- predict(gfit, n.ahead=1)[1]
  df <- coef(gfit)["shape"]
  ES <- sigma*(dt(qt(p,df),df)/(1-p))*((df + (qt(p,df))^2)/(df-1))
  return(ES)
}
ESgarch

## Date vectors for backtest

from <-time(SRLOSS)[-c((nrow(SRLOSS)-999): nrow(SRLOSS))]
#from
to <- time( SRLOSS)[-c(1:1000)]
#to
#SREES <- fapply(SRLOSS, from, to, FUN = ESgarch)

SREES <- applySeries(SRLOSS, from=from, to=to, FUN= ESgarch)
SREES 
SRESSL1 <- lag(SREES, k=1)
SRESSL1
res <- na.omit(cbind(SREES,SRESSL1))
res
colnames(res) <- c("SELOSS","ES99")
plot(res[,2], col= "red", ylim= range(res), main= "BK: t-GARCH(1,1) ES 99 %", ylab= "Percentages", xlab= "Time")
points(res[,1],type = "p", cex=0.2, pch= 19, col= "blue")
legend("topleft", legend= c("Loss","ES"), col = c( " blue" , "red" ) , lty = c(NA, 1) , pch = c(19, NA))


##################################################################################################################

### Computing daily return of Bank of Kigali

##################################################################################################################

log_returns <- diff(log(NewT), lag=1)
log_returns
#newtime <- data1[,1]-1
#log_returns1 <- data.frame(data1[,1],log_returns)

plot(log_returns, type="l", col="blue", main="Returns of Bank of Kigali" )
plot(100*log_returns,type="l", xlab="Time",ylab="%",col=1)

 skewness(log_returns)
 kurtosis(log_returns)
 mean(log_returns)
 

 
 
 
### Spectrum
#raw = spectrum(log_returns)
#smooth = spectrum(log_returns,spans=c(25,5,25),main="Smoothed periodogram",ylim=c(1e-5,2e-4))


#histLR <- hist(log_returns, breaks = 70, plot = TRUE, main = "Log Return Distribution")

HistLR <- hist(log_returns, breaks = 50, plot = FALSE)

hchart(HistLR) %>% 
  hc_title(text = "Log Return Distribution")

 hist(log_returns, breaks = 50, plot = TRUE, col="blue")
##################################################################################

## garchfit for for GARCH(1,1) Parameters

##################################################################################
fit = garchFit(formula = ~ garch(1,1), data = NewT)
print(fit)
formula(fit)

## coefficients of GARCH (1,1)
coef(fit)

## Plot of which
#plot(fit)
plot(fit, which = 2)



####################################################################################

###### MODELLING THE VOLATILITY OF STOCK RETURN

## garchfit for for GARCH(1,1) Parameters
#####################################################################################


fitReturn = garchFit(formula = ~ garch(1,1), data = log_returns)
print(fitReturn )
formula(fitReturn )

## Coefficients of GARCH (1,1)
coef(fitReturn )

## Plot of which
#plot(fit)
plot(fitReturn , which = 2)




### Volatility Modelling 
### Standard  Deviation

Volatility1 = volatility(fit, type = "sigma")
head(Volatility1)
class(Volatility1)
View(Volatility1)

### Variance

Volatility2 = volatility(fit, type = "h")
head(Volatility2)
class(Volatility2)
View(Volatility2)

#################################################################################
### IMPORTANT
### VOLATILITY OF STOCK RETURN
#################################################################################


VolatilityR = volatility(log_returns, type = "sigma")
head(VolatilityR)
class(VolatilityR)
View(VolatilityR)
summary(VolatilityR)


### PLOTTING VOLATILIITY of stock returns

plot(VolatilityR, type="l", col="red", main = "Volatility of stock return ", xlab = "Time(Daily)", ylab= "Volatility")



######################################################################
#### MODELLING THE VOLATILIY OF STOCK RETURN WITH GARCH(1,1)

## MY WORK
######################################################################
### Standard  Deviation

VolatilityReturn = volatility(fitReturn, type = "sigma")
head(VolatilityReturn)
class(VolatilityReturn)
View(VolatilityReturn)

### Variance

VolatilityReturn2 = volatility(fitReturn, type = "h")
head(VolatilityReturn2)
class(VolatilityReturn2)
View(VolatilityReturn2)





### PLOTTING VOLATILIITY

plot(VolatilityReturn, type="l", col="red", main = "Volatility of stock return using GARCH(1,1)", xlab = "Time(Daily)", ylab= "Volatility")



######################################################################################################

#### ALTERNATIVE METHOD

####################################################################################################

garch11 = garch(log_returns, order = c(1, 1), series = NULL, control = garch.control())

garch.control(maxiter = 200, trace = TRUE, start = NULL,
              grad = c("analytical","numerical"), abstol = max(1e-20, .Machine$double.eps^2),
              reltol = max(1e-10, .Machine$double.eps^(2/3)), xtol = sqrt(.Machine$double.eps),
              falsetol = 1e2 * .Machine$double.eps)



summary(garch11)
coef(garch11)

plot(garch11)

Volatility11 = volatility(garch11)
head(Volatility11)
class(Volatility11)
View(Volatility11)


library(AER)
data("NYSESW")
head(NYSESW, 10)



arch.test(log_returns, output = TRUE)
jarque.bera.test(log_returns)



x=log_returns
ArchTest <- function (x, lags=12, demean = FALSE) 
{
  # Capture name of x for documentation in the output  
  xName <- deparse(substitute(x))
  # 
  x <- as.vector(x)
  if(demean) x <- scale(x, center = TRUE, scale = FALSE)
  #  
  lags <- lags + 1
  mat <- embed(x^2, lags)
  arch.lm <- summary(lm(mat[, 1] ~ mat[, -1]))
  STATISTIC <- arch.lm$r.squared * length(resid(arch.lm))
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- lags - 1
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "ARCH LM-test;  Null hypothesis:  no ARCH effects"
  result <- list(statistic = STATISTIC, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name =
                   xName)
  class(result) <- "htest"
  return(result)
}

