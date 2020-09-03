
library(network)


setwd("C:/Users/dominik_deffner/Documents/GitHub/CulturalEffectivePopulationSize")

N <- 100
Density <- 0.2

x <- matrix(NA, nrow = N, ncol = N)

for (i in 1:N) {
  for (j in 1:N) {
    if ( i == j){
      x[i,j] <- 1
      x[j,i] <- 1
    } else {
    Tie <- rbinom(1, 1, Density)
    x[i,j] <- Tie
    x[j,i] <- Tie
    }
  }
}




network <- as.matrix(x)





routes_network <- network(network, ignore.eval = FALSE, directed = FALSE, matrix.type = "adjacency")

plot(routes_network, vertex.cex = 3, col="green")

a  <- summary(routes_network)
