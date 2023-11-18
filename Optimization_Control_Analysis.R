############: Optimization control analysis using Pontryagin’s Maximum Principle

# Required libraries
library(deSolve)
library(rootSolve)
library(optimx)

# Defining the system of ordinary differential equations
#Step 1: Construct the model

model3 <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    #model functions
    egg_prod_fn_pc <- function(M,k,gma) {1/(1+(M/k)*(1-exp(-gma)))^(k+1)}
    sex_reprod_fn <- function(M,k,gma) {1-((1+(M/k)*(1-exp(-gma)))/(1+(M/k)*(2-exp(-gma))))^(k+1)}
    egg_prod_fn <- function(M,k,gma) {M * egg_prod_fn_pc(M,k,gma)}
    egg_prod_fn_sr <- function(M,k,gma) {egg_prod_fn(M,k,gma)*sex_reprod_fn(M,k,gma)}
    
    cp <- function(gp,h,tau) {(-log(1-gp*h))/tau}
    cc <- function(gc,h,tau) {(-log(1-gc*h))/tau}
    ca <- function(ga,h,tau) {(-log(1-ga*h))/tau}
    
    #model ODE equations
    dMp = betap*(1-phi)*L - (mu+cp(gp,h,tau))*Mp
    dMc = betac*(1-phi)*L - (mu+cc(gc,h,tau))*Mc
    dMa = betaa*(1-phi)*L - (mu+ca(ga,h,tau))*Ma
    dL  = ((egg_prod_fn_sr(Mp,k,gma)*np*lambdap*(1-phi)) + 
             (egg_prod_fn_sr(Mc,k,gma)*nc*lambdac*(1-phi)) + 
             (egg_prod_fn_sr(Ma,k,gma)*na*lambdaa*(1-phi))) - 
      (muL + ((1-phi)*(betap*np + betac*nc + betaa*na)))*L;
    
    list(c(dMp, dMc, dMa, dL))
  })
}
