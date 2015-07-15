library(minpack.lm)
library(deSolve)

ABCD_First_Order <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dA <- -k1*A
    dB <- k1*A-k2*B
    dC <- k2*B-k3*C
    dD <- k3*C
    list(c(dA, dB, dC, dD))
  })
}

func_ode <- function(params,tt){
  parameters <- c(k1=params$k1, k2=params$k2, k3=params$k3)
  state <- c(A=1,B=0,C=0,D=0)
  out <- ode(y = state, times = tt, func = ABCD_First_Order, parms = parameters)
  fout <- rowSums(out[,3:5])
  return(fout[1:length(tt)])
}

t <- seq(0, 100, by = 0.01)

## values over which to simulate data
pp = list(k1=0.05,k2=0.01,k3=0.01)
simDNoisy <- func_ode(pp,t) + rnorm(10001,sd=.001)

## residual function
residFun <- function(p, observed, tt) observed - func_ode(p,tt)

## starting values for parameters
parStart <- list(k1=0.08,k2=0.01,k3=0.01)
lower <- c(k1=0.0,k2=0.0,k3=0.0)
upper <- c(k1=1,k2=1,k3=1)

## perform fit
nls.out <- nls.lm(par=parStart,lower=lower,upper=upper, fn = residFun, observed = simDNoisy, tt = t, control = nls.lm.control(nprint=1))

## plot model evaluated at final parameter estimates
#lines(x,getPred(as.list(coef(nls.out)), x), col=2, lwd=2)
## summary information on parameter estimates
summary(nls.out)
