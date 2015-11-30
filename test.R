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

parameters <- c(k1 = 0.05, k2 = 0.10, k3 = 0.01)
state <- c(A=1, B=0, C=0, D=0)
times <- seq(0, 100, by = 0.01)
out <- ode(y = state, times = times, func = ABCD_First_Order, parms = parameters)
plot(out)
