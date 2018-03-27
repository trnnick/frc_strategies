
#### ---- Load required packages ---- ####
if (!require("smooth")){install.packages("smooth")};library(smooth)
if (!require("vars")){install.packages("vars")};library(vars)
if (!require("thief")){install.packages("thief")};library(thief)
if (!require("MAPA")){install.packages("MAPA")};library(MAPA)

#### ---- Helper functions ----  ####
# Functions for iterative, direct and pseudo-MIMO
# This section is only useful to get a view on the inner working of the methods, but is not essential,
# as we produce all forecasts using the smooth package as well, that is more general and can handle 
# any ARIMA.

# Function to generate out-of-sample forecasts
ar1.forecast <- function(y,a,h){
  # y is the vector of historical data
  # a is the coefficient of the model
  # h is the forecast horizon
  f <- vector("numeric",h)
  for (i in 1:h){
    if (i == 1){
      f[i] <- tail(y,1)*a
    } else {
      f[i] <- f[i-1]*a
    }
  }
  return(f)
}

# Function to generate in-sample forecasts. This is needed for numerical optimisation.
ar1.insample <- function(y,a,h){
  # y is the vector of historical data
  # a is the coefficient of the model
  # h is the forecast horizon
  n <- length(y)
  f <- array(NA,c(n,h),dimnames=list(NULL,paste0("t+",1:h)))
  for (i in 1:h){
    if (i == 1){
      f[(1+i):n,i] <- a*y[i:(n-1)]
    } else {
      f[(1+i):n,i] <- a*f[i:(n-1),i-1]
    }
  }
  return(f)
}

# Function to optimise AR(1) coefficient.
ar1.optimise <- function(y,h){
  # y is the vector of historical data
  # h is the forecast horizon
  a <- optim(0,ar1.cost,method="Brent",lower=-1,upper=1,y=y,h=h)$par
  return(a)
}

# Cost function, MSE with horizon argument.
ar1.cost <- function(a,y,h){
  # a is the coefficient of the model
  # y is the vector of historical data
  # h is the forecast horizon
  f <- ar1.insample(y,a,max(h))
  f <- f[,h,drop=FALSE]
  se <- sum(colMeans((matrix(rep(y,length(h)),ncol=length(h)) - f)^2,na.rm=TRUE))
  return(se)
}

# Analytical parameter estimation.
ar1.estimate <- function(y,h){
  n <- length(y)
  yl <- matrix(y[1:(n-h)],ncol=1)
  yy <- matrix(y[(h+1):n],ncol=1)
  a <- solve(t(yl) %*% yl) %*% t(yl) %*% yy
  return(a)
}

#### ---- Demonstrate forecasting strategies ---- ####
# Here we demonstarte the different strategies
# We track estimated parameters and forecasts

# For this example we will use a forecast horizon of 3
h <- 3

# Simulate AR(1) using sim.ssarima() from smooth
set.seed(123456)
n <- 36           # Time series sample size
# We also remove firts 12 values, as a burn-in period
y <- sim.ssarima(orders=list(ar=1), lags=1, obs=n+12, frequency=12, AR=0.5)$data[-(1:12)]

# Separate time series into in-sample and test set
y.in <- y[1:(n-h)]
y.out <- y[(n-h+1):n]

# Iterative approach
# Get parameters
a1 <- ar1.optimise(y.in,1)          # Get parameter by numerical optimisation
a1e <- ar1.estimate(y.in,1)         # or analytically. 
a1 - a1e                            # There are minor differences as expected.

# Forecast
frc1 <- ar1.forecast(y.in,a1,h)
frc1e <- ar1.forecast(y.in,a1e,h)
# Or alternatively we use the (state-space) ARIMA implementation in the smooth
# package that allows various cost functions
model1 <- ssarima(y.in,orders=list(ar=1), lags=1, cfType="MSE")
a1m <- as.vector(model1$persistence)
frc1m <- forecast(model1,h=h)$mean

#### Direct approach ####
# A similar logic is followed, but now h models are created.
a2 <- a2e <- frc2 <- frc2e <- vector("numeric",h)
for (i in 1:h){
  a2[i] <- ar1.optimise(y.in,i)    # These parameters are for producing t+h forecasts, i.e. a normal AR(1) model optimal for t+h forecasts.
  a2e[i] <- ar1.estimate(y.in,i)   # These parameters are for producing t+1 forecasts directly
  frc2[i] <- ar1.forecast(y.in,a2[i],i)[i]  # t+h direct forecast equivalents
  frc2e[i] <- ar1.forecast(y.in,a2e[i],1)   # t+1 direct forecasts
}
# To compare the parameters discount a2 for the h-steps ahead
a2^(1:3) - a2e

# The differences are small, as expected, but qualitatively the two approaches are the same.
# Using the smooth package we need to set cfType="MSEh" and h for the appropriate horizon.
model2 <- model2e <- list()
a2m <- a2me <- frc2m <- vector("numeric",h)
for (i in 1:h){
  model2[[i]] <- ssarima(y.in,orders=list(ar=1), lags=1, cfType="MSEh", h=i)
  a2m[i] <- as.vector(model2[[i]]$persistence) # This is the equivalent to a2, but different values do to differences in implementation, initialisation
  frc2m[i] <- forecast(model2[[i]],h=i)$mean[i]
  # Or we can implement the same with conventional t+1 cost, but increasing the lag order
  model2e[[i]] <- ssarima(y.in,orders=list(ar=1), lags=i, cfType="MSE")
  a2me[i] <- as.vector(model2e[[i]]$persistence)[i] # This is the equivalent to a2e.
}

#### Pseudo-MIMO ####
a3 <- ar1.optimise(y.in,1:h)
# There is no analytical equiavlent and now we rely only on optimising across all h-steps
frc3 <- ar1.forecast(y.in,a3,h)
# Alternatively, using the smooth package we need to set cfType="TMSE"
model3 <- ssarima(y.in,orders=list(ar=1), lags=1, cfType="TMSE",h=h)
a3m <- as.vector(model3$persistence)
frc3m <- forecast(model3,h=h)$mean

#### MIMO ####
# Instead of implementing VAR(1) from scratch we will rely on a typical VAR implementation
# First we create the multivariate time series, including lags up to the horizon-1
Z <- array(NA,c(length(y),h),dimnames=list(NULL,c("Y",paste("Ylag",1:(h-1)))))
for (i in 1:h){
  Z[i:n,i] <- y[1:(n-i+1)]
}
# Remove h-1 first rows that contain NAs
Z <- Z[h:n,]
Z <- ts(Z,frequency=frequency(y),end=end(y))
# Estimate the unrestricted model - This is not what we want, as each row includes inputs
# from the lagged variables as well
varModel <- VAR(Z, p=1, type="none")
# Create restriction
rest <- array(0,c(h,h))
rest[,h] <- 1
# Re-estimate the restricted model - This is MIMO forecast
varModelRest <- restrict(varModel, method="man", resmat=rest)
print(varModelRest)
# The output should be read in the following way
# Estimated coefficients for equation Y --> That is the t+h forecast
# Estimated coefficients for equation Ylag.1 --> That is the t+h-1 forecast
# ...
# Get model parameters
a4 <- as.vector(unlist(lapply(varModelRest$varresult,
                              function(x){x$coefficients}))[h:1]) # These should be compared
                                                # with a2e and a2me, which ignore covariances
# Get forecasts
varForecast <- forecast(varModelRest,h=h)
frc4 <- lapply(varForecast$forecast,function(x){as.vector(x$mean)})[h:1] # Re-order from t+1 to t+h
frc4 <- as.vector(diag(matrix(unlist(frc4),nrow=h,byrow=TRUE)))

#### Measure accuracy ####
MSE <- colMeans((matrix(rep(y.out,4),ncol=4) - cbind(frc1e,frc2e,frc3,frc4))^2)
names(MSE) <- c("Iterative","Direct","Pseudo-MIMO","MIMO")
round(MSE,4)

# MSE for the forecasts produced using smooth functions
MSEsmooth <- colMeans((matrix(rep(y.out,4),ncol=4) - cbind(frc1m,frc2m,frc3m,frc4))^2)
names(MSEsmooth) <- c("Iterative","Direct","Pseudo-MIMO","MIMO")
round(MSEsmooth,4)

# Finally, we conclude the note by making the point that temporal aggregation has some connections 
# with these forecasting strategies. For comparison, here are the MAPA and THieF forecasts
y.in <- ts(y.in,frequency=12,start=c(1,1))
MSEaggr <- c(mean((y.out - thief(y.in,h=3)$mean[1:3])^2),
             mean((y.out - mapa(y.in,fh=3)$outfor)^2))
names(MSEaggr) <- c("THieF","MAPA")
round(MSEaggr,4)
