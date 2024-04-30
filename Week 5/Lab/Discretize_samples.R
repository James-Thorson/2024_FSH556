
setwd( R'(C:\Users\James.Thorson\Desktop\Git\2024_FSH556_private\Week 5\Lab)' )
vismba = readRDS( "vismba.rds" )

#       
library(sf)
library(fmesher)
samples = data.frame("x"=vismba$gx, "y"=vismba$gy, "agb"=vismba$agb )
samples = st_as_sf( samples, coords=c("x","y") )
#samples = samples[1:2,,drop=FALSE]
samples$agb = 1

grid = st_make_grid( st_bbox(c(xmin=0, xmax=1000, ymin=0, ymax=500)), n=c(8,4) )
grid_i = st_intersects( samples, grid )
Count_i = tapply( samples$agb, INDEX=factor(unlist(grid_i),levels=1:length(grid)), FUN=length )
Data = data.frame( st_coordinates(st_centroid(grid)), "Count"=ifelse(is.na(Count_i),0,Count_i) )

grid_sf = st_sf(grid, Count=Count_i)
plot( grid_sf, axes=TRUE, reset=FALSE, pal=sf.colors(n=10, alpha=0.2), breaks=seq(0,max(Data$Count),length=11) )
plot( samples, add=TRUE, pch=20 )

# Get adjacency
st_rook = function(m, ...) st_relate(m, m, pattern="F***1****", ... )
grid_A = st_rook( grid_sf, sparse=TRUE )
A = as(grid_A,"sparseMatrix")

# Get SPDE matrices
mesh = fm_mesh_2d( st_coordinates(st_centroid(grid_sf)), 
                   plot.delay = NULL, 
                   refine = TRUE )
spde <- fm_fem( mesh, 
                order = 2 )

# M0 = spde$c0
# M1 = spde$g1
# M2 = spde$g2
