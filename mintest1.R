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

ABC_First_Order <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dA <- -k1*A
    dB <- k1*A-k2*B
    dC <- k2*B
    list(c(dA, dB, dC))
  })
}

AB_First_Order <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dA <- -k1*A
    dB <- k1*A
    list(c(dA, dB))
  })
}

func_ode_2 <- function(params,tt){
  parameters <- c(k1=params$k1)
  state <- c(A=1,B=0)
  out <- ode(y = state, times = tt, func = AB_First_Order, parms = parameters)
  fout <- out[,3]
  return(fout[1:length(tt)])
}

func_ode_3 <- function(params,tt){
  parameters <- c(k1=params$k1, k2=params$k2)
  state <- c(A=1,B=0,C=0)
  out <- ode(y = state, times = tt, func = ABC_First_Order, parms = parameters,method="ode45")
  fout <- rowSums(out[,3:4])
  return(fout[1:length(tt)])
}

func_ode_4 <- function(params,tt){
  parameters <- c(k1=params$k1, k2=params$k2, k3=params$k3)
  state <- c(A=1,B=0,C=0,D=0)
  out <- ode(y = state, times = tt, func = ABCD_First_Order, parms = parameters)
  #fout <- rowSums(out[,3:5])
  #return(fout[1:length(tt)])
  return(out[,3:5])
}


t <- seq(0, 50, by = 0.1)

## values over which to simulate data
pp = list(k1=0.99,k2=0.317,k3=0.08)
simDNoisy <- func_ode_4(pp,t) #+ rnorm(length(t),sd=.01)
dat <- data.frame(simDNoisy)

plot(x=t,dat$D)


## residual function
residFun <- function(p, observed, tt) {
  residFun <- func_ode_3(p,tt) - observed
}
## starting values for parameters
parStart <- list(k1=0.6,k2=-0.5)
lower <- c(k1=0.5,k2=-0.4)
upper <- c(k1=0.9,k2=-0.01)

## perform fit
nls.out <- nls.lm(par=parStart,lower=lower, upper=upper, fn = residFun, observed = simDNoisy, tt = t, control = nls.lm.control(nprint=1))

## plot model evaluated at final parameter estimates
#plot(t,getPred(as.list(coef(nls.out)), x), col=2, lwd=2)
layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,1.4), respect = FALSE)
par(mar = c(0, 4.1, 3, 2.1))
plot(t,simDNoisy)
lines(t,func_ode_2(nls.out$par,t),col="red",xaxt="n")
par(mar = c(4.1, 4.1, 0.3, 2.1))
plot(t,nls.out$fvec)

## summary information on parameter estimates
summary(nls.out)
