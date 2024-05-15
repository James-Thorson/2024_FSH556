
data_dir = "C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Chap_8/"

setwd( R'(C:\Users\James.Thorson\Desktop\Git\2024_FSH556_private\Week 8\Lab)' )

# Libraries and functions
library( TMB )
library( fmesher )
library( sf )
source("C:/Users/James.Thorson/Desktop/Git/Spatio-temporal-models-for-ecologists/Shared_functions/add_legend.R")

# 
if( FALSE ){
  EBS = st_read( "C:/Users/James.Thorson/Desktop/Git/FishStatsUtils/inst/region_shapefiles/EBSshelf/EBSshelf.shp" )
  saveRDS( EBS, file.path( R'(C:\Users\James.Thorson\Desktop\Git\2024_FSH556_private\Week 8\Lab)', "EBS.rds") )
}

#
pollock = readRDS( file.path(data_dir,"pollock.rds") )
pollock = st_as_sf(pollock, coords = c("Long","Lat"), crs="+proj=longlat +datum=WGS84" )
pollock = st_transform( pollock, crs=st_crs("+proj=utm +zone=2 +datum=WGS84 +units=km") )

#
survey_domain = readRDS( "EBS.rds" )
survey_domain = st_transform( survey_domain, crs=st_crs(pollock) )

# Reduce data to EBS domain
pollock = st_intersection( pollock, survey_domain )

# Make triangulated mesh
mesh = fm_mesh_2d( st_coordinates(pollock), cutoff=100, refine=TRUE )
# Create matrices in INLA
spde = fm_fem(mesh, order=2)
# create projection matrix from vertices to samples
A_is = fm_evaluator( mesh, loc=st_coordinates(pollock) )$proj$A
# Create extrapolation grid
cellsize = 25
grid = st_make_grid( survey_domain, cellsize=cellsize )
grid = st_intersection( grid, survey_domain )
# create projection matrix from vertices to grid
A_gs = fm_evaluator( mesh, loc=st_coordinates(st_centroid(grid)) )$proj$A

# Compile
Version = "pollock_index"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

# Make inputs
year_set = min(pollock$Year):max(pollock$Year)
Data = list( "n_t" = length(year_set),
             "a_g" = as.numeric(st_area(grid)),
             "z_g" = st_coordinates(st_centroid(grid))[,2],
             "b_i" = pollock$Wt,
             "a_i" = pollock$AreaSwept_ha / 100, # Convert hectares to km^2
             "t_i" = pollock$Year - min(pollock$Year),
             "A_is" = A_is,
             "A_gs" = A_gs,
             "M0" = spde$c0,
             "M1" = spde$g1,
             "M2" = spde$g2 )
Params = list( "beta_t" = rep(0,Data$n_t),
               "ln_tauO" = log(1),
               "ln_tauE" = log(1),
               "ln_kappa" = 1,
               "ln_phi" = log(1),
               "logit_rhoE" = 0,
               "finv_power" = 0,
               "omega_s" = rep(0, mesh$n),
               "epsilon_st" = matrix(0, nrow=mesh$n, ncol=Data$n_t) )
Random = c("omega_s", "epsilon_st")

# Build and run
Obj = MakeADFun( data=Data, parameters=Params, random=Random )
Obj$env$beSilent()

# Run
Opt = nlminb( start=Obj$par, obj=Obj$fn, gr=Obj$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4) )
Report = Obj$report()
Opt$SD = sdreport( Obj, 
                   bias.correct=FALSE, 
                   getJointPrecision=TRUE )



