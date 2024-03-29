
setwd( R'(C:\Users\James.Thorson\Desktop\Git\2024_FSH556\Week 1\Lab)' )

mu = 2
n_indiv = 30
#n_obsperindiv = 3
#sd_indiv = 1

#eps_j = rnorm( n_indiv, mean=mu, sd=sd_indiv )

Y_i = rpois( n_indiv, mu )

Data = list( Y_i=Y_i )
Params = list( "log_mu"=0 )

library(TMB)
compile("poisson_glm.cpp")
dyn.load( dynlib("poisson_glm") )

obj = MakeADFun( parameters=Params, 
                 data = Data )
opt = nlminb( start=obj$par, 
              obj = obj$fn, 
              gr = obj$gr )

Lm = glm( Y_i ~ 1, family=poisson(link="log") )
