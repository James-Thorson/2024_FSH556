
setwd( R'(C:\Users\James.Thorson\Desktop\Git\2024_FSH556\Week 1\Homework)' )

library(VAST)
data( EBS_pollock_data )
CPUE = EBS_pollock_data[[1]]
write.csv( CPUE, "EBS_pollock_data.csv", row.names=FALSE )
