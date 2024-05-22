
###########################
# One dimensional
###########################

#library(Matrix)
n_g = 200

# Domain characteristics
lat_g = seq(55, 65, length=n_g)
Temp_g = seq(15, 5, length=n_g)

# Parameters
# if testing diffusion via comparison with variance of displacement
unittime_variance = 2^2
D = unittime_variance / 2
Dprime = D / mean(diff(lat_g))^2
preference_g = -0 * (Temp_g - 10)^2
# if testing dynamics for taxis
  #diffusion_coefficient = 1 ^ 2
  #preference_g = -0.1 * (Temp_g - 10)^2

# Movement operator
A_gg = ifelse( round(abs(outer(lat_g,lat_g,"-")),2) == round(mean(diff(lat_g)),2), 1, 0 )
# Diffusion
diffusion_gg = A_gg * Dprime
diag(diffusion_gg) = -1 * colSums(diffusion_gg)
# Taxis
taxis_gg = A_gg * outer(preference_g, preference_g, "-")
diag(taxis_gg) = -1 * colSums(taxis_gg)
# Total
mrate_gg = diffusion_gg + taxis_gg
if( any(mrate_gg-diag(diag(mrate_gg))<0) ) stop("`mrate_gg` must be a Metzler matrix. Consider changing parameterization")
# Annualized
mfraction_gg = Matrix::expm(mrate_gg)

# plot
matplot( scale(cbind(lat_g,Temp_g,preference_g,mfraction_gg[,ceiling(0.25*n_g)],mfraction_gg[,ceiling(0.75*n_g)])), type="l", lty="solid", lwd=2 )

#
Matrix::image(mfraction_gg, xlab="From", ylab="To", border=NA)

# Stationary distribution
stationary_g = eigen(mfraction_gg)$vectors[,1]
stationary_g = stationary_g / sum(stationary_g)

#
matplot( x=lat_g, y=cbind(preference_g-min(preference_g),stationary_g), type="l", col=c("black","blue") )

# n(t+1) = Mrate_gg * n(t)

# Diffusion should equal variance of displacement
  # Using central cell to minimize boundary issues
# sum( mfraction_gg[,6] * 1:n_g )
mean1 = weighted.mean(lat_g, w=mfraction_gg[,ceiling(n_g/2)])
# sum( mfraction_gg[,6] * (1:n_g-mean1)^2 )
var1 = weighted.mean( (lat_g-mean1)^2, w=mfraction_gg[,ceiling(n_g/2)])
# Should be equal
var1
2 * 1 * D  # MSD = 2nDt

#######################
# 2-dimensional
#######################

library(sf)

# Define parameters
bounds = 20
shape_type = c("square", "hex")[1]
cellsize = 0.8
n_neighbors = ifelse(shape_type=="square", 4, 6)
n_dim = 2  # 2-dimensions
unittime_variance = 1^2

# Define dimain and discretization
sf_area = st_polygon(x = list(cbind(c(0,0,bounds,bounds,0),c(0,bounds,bounds,0,0))) )
sf_grid = st_make_grid( sf_area, cellsize=cellsize, square = isTRUE(shape_type=="square") )

# Get adjacency
st_rook = function(a, b = a, ...) st_relate(a, b, pattern = "F***1****", ... )
grid_A = st_rook( sf_grid, sparse=TRUE )
A_gg = as(grid_A,"sparseMatrix")

# Calculcate D
D = unittime_variance / 2
# Dprime presumably depends on discretization size and shape
#Dprime = unittime_variance / n_neighbors / cellsize^2
Dprime = D / cellsize^2

# Assemble diffusion-rate matrix
diffusion_gg = A_gg * Dprime
diag(diffusion_gg) = -1 * Matrix::colSums(diffusion_gg)
# integrated movement fractions matrix
mfraction_gg = Matrix::expm(diffusion_gg)

# Using central cell to minimize boundary issues
sf_point = st_point(c(10.1,10.1) )
centroid = as.numeric( st_within(sf_point,sf_grid) )
# Calculate average in X and Y dimensions
meanX = weighted.mean( x=st_coordinates(st_centroid(sf_grid))[,1], w=mfraction_gg[,centroid])
meanY = weighted.mean( x=st_coordinates(st_centroid(sf_grid))[,2], w=mfraction_gg[,centroid])
# calculate variance in X and Y dimensions
varX = weighted.mean( (st_coordinates(st_centroid(sf_grid))[,1]-meanX)^2, w=mfraction_gg[,centroid])
varY = weighted.mean( (st_coordinates(st_centroid(sf_grid))[,2]-meanY)^2, w=mfraction_gg[,centroid])
# Mean-squared-displacement (MSD) = var1 + var2
(MSD = varX + varY)
#2 * D   # MSD = 2nD from https://en.wikipedia.org/wiki/Mean_squared_displacement
2 * n_dim * D    # Only works for square cells

##################
# Analytical check
##################

## Fundamental solution
set.seed(1)
d <- function(x,t) mvtnorm::dmvnorm(x, mean=rep(0,n), sigma=Sigma*t)

## Spatial dimension
n <- 3

## Random covariance matrix
Chol = matrix(rnorm(n*n), nrow=n, ncol=n)
Sigma <- Chol %*% t(Chol)

## Random spatial coordinate
x0 <- abs(rnorm(n))

## Random time
t0 <- abs(rnorm(1))

## Check diffusion equation
numDeriv::grad(function(t) d(x0, t), t0) ## LHS (time derivative)
.5*sum(Sigma*numDeriv::hessian(function(x)d(x, t0), x0)) ## RHS (space derivative)

