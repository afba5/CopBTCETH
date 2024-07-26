library(quantmod)
library(qrmtools)

df_BTC <- getSymbols('BTC-EUR',src='yahoo',auto.assign=FALSE, from = "2017-11-11",
                     to = "2024-05-10")
df_ETH <- getSymbols('ETH-EUR',src='yahoo',auto.assign=FALSE, from = "2017-11-11",
                     to = "2024-05-10")

BTC_ETH <- cbind(df_BTC$`BTC-EUR.Adjusted`,df_ETH$`ETH-EUR.Adjusted`)
colnames(BTC_ETH) <- c("BTC-EUR","ETH-EUR") 

R_BTCETH <- returns(BTC_ETH) # calculate log-returns 

library(rugarch)

# For BTC:
# GARCH Model	: eGARCH(3,4) 
# Mean Model	: ARFIMA(0,0,0) WITH MEAN
# Distribution	: sstd

ctrl = list(RHO = 1,DELTA = 1e-8,MAJIT = 100,MINIT = 650,TOL = 1e-6)
## Specify model: ARMA(0,0)-eGARCH(3,4) WITH MEAN and skew Student t innovations
spec_BTC = ugarchspec(variance.model = list(model="eGARCH",garchOrder = c(3, 4)),mean.model = list(armaOrder = c(0,0),include.mean=TRUE),distribution.model = "sstd") 
## Fit marginal ARMA-GARCH model for BTC-EUR returns
garch.fit_BTC = ugarchfit(data = R_BTCETH$`BTC-EUR`, spec = spec_BTC, solver = "solnp", solver.control = ctrl)
## Extract standardized residuals
residuals_BTC <- residuals(garch.fit_BTC,standardize=TRUE)

# For ETH:
# GARCH Model	: eGARCH(4,3) 
# Mean Model	: ARFIMA(0,0,0) WITH ZERO MEAN
# Distribution	: sstd 

## Specify model: ARMA(0,0)-eGARCH(4,3) WITH ZERO MEAN and skew Student t innovations
spec_ETH = ugarchspec(variance.model = list(model="eGARCH",garchOrder = c(4, 3)),mean.model = list(armaOrder = c(0,0),include.mean=FALSE),distribution.model = "sstd") 
## Fit marginal ARMA-GARCH model for ETH-EUR returns
garch.fit_ETH = ugarchfit(data = R_BTCETH$`ETH-EUR`, spec = spec_ETH, solver = "solnp", solver.control = ctrl)
## Extract standardized residuals
residuals_ETH <- residuals(garch.fit_ETH,standardize=TRUE)

# Put together stand. residuals from BTC and ETH ARMA-GARCH models
residuals_BTCETH <- cbind(residuals_BTC,residuals_ETH)
colnames(residuals_BTCETH) <- c("residuals_BTC","residuals_ETH")
# Convert to data frame
df_residuals_BTCETH <- fortify.zoo(residuals_BTCETH) 

library(copula)

# Transform standardized residuals into observations in the interval [0,1]
U <- pobs(df_residuals_BTCETH[,-1]) 

# Time-invariant Gumbel Survival Copula or rotated Gumbel Copula
fitcop_gumbel_s <- fitCopula(rotCopula(gumbelCopula(dim = 2)), data = U, method = "mpl")

# Time-varying rotated Gumbel Copula (180 degrees) or Gumbel Survival

## Calculate the negative copula log-likelihood of Gumbel copula with time-varying parameter 
gumbel_tvp1_CL <- function(theta, data, kappabar) {
  
  ### This function is a modified version in R from Andrew J. Patton's Matlab code
  ### Date: Sunday, 5 Aug, 2001.
  ### Source: https://public.econ.duke.edu/~ap172/code.html
  
  T <- nrow(data)
  u <- data[, 1]
  v <- data[, 2]
  
  kappa <- rep(-999.99, T)
  kappa[1] <- kappabar
  for (jj in 2:T) {
    if (jj <= 10) {
      psi1 <- theta[1] + theta[2] * kappa[jj - 1] + theta[3] * mean(abs(u[1:(jj - 1)] - v[1:(jj - 1)]))
    } else {
      psi1 <- theta[1] + theta[2] * kappa[jj - 1] + theta[3] * mean(abs(u[(jj - 10):(jj - 1)] - v[(jj - 10):(jj - 1)]))
    }
    kappa[jj] <- 1.0001 + psi1^2
  }
  
  ut <- -log(data[, 1])
  vt <- -log(data[, 2])
  
  CL <- -(ut^kappa + vt^kappa)^(1 / kappa) - log(u) - log(v)
  CL <- CL + (kappa - 1) * (log(ut) + log(vt)) - (2 - 1 / kappa) * (log(ut^kappa + vt^kappa))
  CL <- CL + log((ut^kappa + vt^kappa)^(1 / kappa) + kappa - 1)
  CL <- sum(CL)
  CL <- -CL
  
  if (!is.double(theta)) {
    CL <- 1e6
  } else if (!is.double(CL)) {
    CL <- 1e7
  } else if (is.nan(CL)) {
    CL <- 1e8
  } else if (is.infinite(CL)) {
    CL <- 1e9
  }
  
  return(list(CL = CL, kappa = kappa))
}

## Create a second function equal to the one before (gumbel_tvp1_CL) to be used for optimization
gumbel_tvp1_CL_only <- function(theta, data, kappabar) { 
  CL <- gumbel_tvp1_CL(theta, data, kappabar)[1]
  return(CL)
}

## Maximize the copula-likelihood function
alpha1_tvgumbel_s <- coef(fitcop_gumbel_s)
theta0 <- c(sqrt(alpha1_tvgumbel_s-1), 0, 0)

result_tvgumbel_s <- optim(
  par = theta0,
  fn = gumbel_tvp1_CL_only,
  control=list(trace=TRUE, maxit=1000),
  data = 1-U,  # we use 1-U since we look for the rotated copula
  kappabar = alpha1_tvgumbel_s
)

## Extract the coefficients that maximize the copula-likelihood function  
theta_tvgumbel_s <- result_tvgumbel_s$par

## Obtain the time-varying copula parameter
alphat_tvgumbel_s <- gumbel_tvp1_CL(theta_tvgumbel_s, 1-U, alpha1_tvgumbel_s)[[2]]

## Plot constant and time-varying copula parameter
library(ggplot2)
T_alphat_tvgumbel_s <- length(alphat_tvgumbel_s)
data_tvgumbel_s <- data.frame(time = df_residuals_BTCETH[,1], 
                              alphat_tvgumbel_s = alphat_tvgumbel_s, 
                              alpha1_tvgumbel_s = alpha1_tvgumbel_s * rep(1,T_alphat_tvgumbel_s))
ggplot(data_tvgumbel_s, aes(x = time)) +
  geom_line(aes(y = alphat_tvgumbel_s, color = "Time-varying")) +
  geom_line(aes(y = alpha1_tvgumbel_s, color = "Constant"), linetype = "dashed") +
  labs(x = "time", y = expression(theta), color = "") +
  theme_classic() +
  scale_color_manual(values = c("red","black")) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.background = element_blank(),
        legend.position = c(0.9,0.97))

