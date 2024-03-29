
setwd( R'(C:\Users\James.Thorson\Desktop\Git\2024_FSH556\Week 1\Lab)' )

mu = 2
n_indiv = 30
n_obsperindiv = 3
sd_indiv = 0.5
j_i = rep( 1:n_indiv, each=n_obsperindiv )

eps_j = rnorm( n_indiv, mean=mu, sd=sd_indiv )
mu_i = rep( eps_j, each=n_obsperindiv )

Y_i = rpois( n_indiv*n_obsperindiv, mu_i )

library(TMB)

# GLM

Data = list( Y_i = Y_i )
Params = list( "log_mu"=0 )

compile("poisson_glm.cpp")
dyn.load( dynlib("poisson_glm") )

obj = MakeADFun( parameters=Params, 
                 data = Data )
opt = nlminb( start=obj$par, 
              obj = obj$fn, 
              gr = obj$gr )

Lm = glm( Y_i ~ 1, family=poisson(link="log") )

# GLMM
Data = list( Y_i = Y_i,
             j_i = j_i - 1 )
Params = list( "log_mu"=0,
               "log_sd_indiv" = 0,
               "eps_j" = rep(0,n_indiv) )

compile("poisson_glmm.cpp")
dyn.load( dynlib("poisson_glmm") )

obj = MakeADFun( parameters=Params, 
                 data = Data,
                 random = "eps_j",
                 DLL = "poisson_glmm" )
opt = nlminb( start=obj$par, 
              obj = obj$fn, 
              gr = obj$gr )
SD = sdreport( obj )

library(lme4)
Lmer = glmer( Y_i ~ (1 | j_i), family=poisson(link="log") )
