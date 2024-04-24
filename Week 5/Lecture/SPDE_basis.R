
library(sf)
library(fmesher)

# Make grid
sf_box = st_polygon( list(cbind(c(0,1,1,0,0),c(0,0,1,1,0))) )
sf_grid = st_make_grid( sf_box, n=c(100,100) )
grid_xy = st_coordinates(st_centroid(sf_grid))

# SPDE mesh
mesh = fm_mesh_2d( grid_xy, 
                   cutoff = 0.20 )

# SPDE basis functions
A_is = fm_evaluator( mesh, 
                     loc = grid_xy )$proj$A

# Plot largest basis (to be easily visible)
which_basis = which.max( colSums(as.matrix(A_is)) )
basis = ifelse( A_is[,which_basis]==0, NA, A_is[,which_basis] )
plotgrid = st_sf( sf_grid, 
                  basis = basis )
plot( fm_as_sfc(mesh) )
plot( plotgrid, 
      border = NA,
      nbreaks = 100,
      add = TRUE )
