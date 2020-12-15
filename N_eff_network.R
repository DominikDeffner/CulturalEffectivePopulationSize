
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission

# We simulate 2 different populations, one starting with all the same variants, one with all unique variants.
# We measure diversity and run simulations until diversities cross indicating they've reached equilibrium, then switch transmission mode


# Evolution on different graphs

# Function to calculate effective population sizes
library(vegan)
library(igraph)
library(ggnetwork)

N_e <- function(k_bar,V_k){
  (N * k_bar - 1 ) / ( (V_k / k_bar) + k_bar  -1 )
}

N = 1000          # Census population size
tmax = 300       # Number of timesteps / generations
mu = 1e-4       # Innovation rate
type = "random" # Network type: "random", "scalefree" or "smallworld"
p = 0.4

# Initialize population with cultural traits

Pop <- matrix(NA, 2, N)

# Maximally diverse
Pop[1,] <- sample(1:N)

#All same
Pop[2,] <- rep(1, N)

# Create numerator for cultural variants in both populations
Counter <- c(max(Pop[1,]), max(Pop[2,]))

# Calculate Simpson diversity in both populations
Div <- c()
for (pop_id in 1:2) {
  Div[pop_id] <- vegan::diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
}

d1 <- Div[1]
d2 <- Div[2]
# Burn-in to reach equilibrium
while(Div[1] > Div[2]){
  Copied <- matrix(NA, 2, N)
  
  for (pop_id in 1:2) {
    #Cultural Transmission
    
    Copied[pop_id, ] <- sample(1:N, N, replace = TRUE)
    Pop[pop_id,] <- Pop[pop_id,Copied[pop_id, ]]
    
    # Innovation
    Innovators <- rbinom(N,1,mu)
    Pop[pop_id, Innovators == 1] <- (Counter[pop_id] + 1) : (Counter[pop_id] + length(which(Innovators==1)))
    Counter[pop_id] <- max(Pop[pop_id,])
  }#pop_id
  
  # Calculate Simpson diversity in both populations
  Div <- c()
  for (pop_id in 1:2) {
    Div[pop_id] <- vegan::diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
  }
  
d1 <- c(d1, Div[1])
d2 <- c(d2, Div[2])

}#while


# Create output objects
N_e <- function(k_bar,V_k){
  (N * k_bar - 1 ) / ( (V_k / k_bar) + k_bar  -1 )
}
#Frequency_spectra <- array(NA,  dim = c(2, tmax, N))
N_effective <- list()

Offspring_Record1 <-  sapply(1:N, function(x) length(which(Copied[1,] == x)))
Offspring_Record2 <-  sapply(1:N, function(x) length(which(Copied[2,] == x)))

N_effective[[1]] <- N_e(mean(Offspring_Record1), var(Offspring_Record1))
N_effective[[2]] <- N_e(mean(Offspring_Record2), var(Offspring_Record2))


# Simulate random notwork
#g <- erdos.renyi.game(N, p, type = "gnp")

# Scale free network
#g <- sample_pa(N, power = 1, m = 1, out.dist = NULL, out.seq = NULL, out.pref = FALSE, zero.appeal = 1, directed = FALSE, algorithm ="psumtree", start.graph = NULL)

# Small world
g <- watts.strogatz.game(1, N, 1, p = 0.05, loops = FALSE, multiple = FALSE)

A <- as.matrix(get.adjacency(g, type = "both"))


for (t in 1:tmax) {
  
  for (pop_id in 1:2) {
    
    
     
    Pop_new <- c()
    Copied <- c()
    
    #Cultural Transmission in network
    for (i in 1:N) {
     # if ( length(which(A[i,] == 1)) >  1 ){
        Chosen <- sample(which(A[i,] == 1), 1 )
     # } else if (length(which(A[i,] == 1)) == 1){
     #   Chosen <- which(A[i,] == 1)
    #  } else {
    #    Chosen <- i
    #  }
      Copied <- c(Copied, Chosen)
      Pop_new[i] <- Pop[pop_id, Chosen]  
    }
    
    # Replace population with juveniles
    Pop[pop_id , ] <- Pop_new

    # Compute effective population size
    Offspring_Record <-  sapply(1:N, function(x) length(which(Copied == x)))
    N_effective[[pop_id]] <-  c(N_effective[[pop_id]], N_e(mean(Offspring_Record), var(Offspring_Record)))
    
    
    # Innovation
    Innovators <- rbinom(N,1,mu)
    Pop[pop_id, Innovators == 1] <- (Counter[pop_id] + 1) : (Counter[pop_id] + length(which(Innovators==1)))
    Counter[pop_id] <- max(Pop[pop_id,])
    

  }#pop_id
  
  print(t)
}




plot(N_effective[[1]], type = "b")



















