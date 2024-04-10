
setwd( R'(C:\Users\James.Thorson\Desktop\Git\2024_FSH556\Week 3\Lab)' )

library(tinyVAST)
data(salmon_returns)

Sockeye = subset( salmon_returns, Species=="sockeye" )
write.csv( Sockeye, "Sockeye_returns.csv", row.names=FALSE )

Y_tc = tapply( Sockeye$Biomass, INDEX=list(Sockeye$Year,Sockeye$Region), FUN=mean )
matplot( log(Y_tc), type="l", lwd=2 )

