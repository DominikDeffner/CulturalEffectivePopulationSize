
library(igraph)
library(ggnetwork)

setwd("C:/Users/dominik_deffner/Documents/GitHub/CulturalEffectivePopulationSize")

N <- 100


# Random networks

p <- 0.01
g <- erdos.renyi.game(N, p, type = "gnp")
A <- as.matrix(get.adjacency(g, type = "both"))
all(sapply(1:N, function(x) length(which(A[x,]==1) )) > 0 )



plot(g, vertex.label= NA, edge.arrow.size=0.5,vertex.size = 1, xlab = "Random network model")


# Scale free networks

g <- sample_pa(N, power = 1, m = 1, out.dist = NULL, out.seq = NULL,
               out.pref = FALSE, zero.appeal = 1, directed = FALSE,
               algorithm ="psumtree", start.graph = NULL)
dens(degree.distribution(g))


plot(g, vertex.label= NA, edge.arrow.size=0.5,vertex.size = 1, xlab = "Scale free network model")


# Small world networks

g <- watts.strogatz.game(1, N, 200, p = 0.01, loops = FALSE, multiple = FALSE)

plot(g, vertex.label= NA, edge.arrow.size=0.5,vertex.size = 1, xlab = "Small world network model")




A <- as.matrix(get.adjacency(g, type = "both"))

sum(A)



