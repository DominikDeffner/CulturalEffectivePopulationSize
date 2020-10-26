
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission

# We simulate 2 different populations, one starting with all the same variants, one with all unique variants.
# We measure diversity and run simulations until diversities cross indicating they've reached equilibrium, then switch transmission mode


#Connectedness

# Function to calculate effective population sizes
library(vegan)

N_e <- function(k_bar,V_k){
  (N * k_bar - 1 ) / ( (V_k / k_bar) + k_bar  -1 )
}

N = 1000          # Census population size
tmax = 300       # Number of timesteps / generations
mu = 1e-3       # Innovation rate
m = 0         #Migration rate between 2 populations
e = 0.3
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
  Div[pop_id] <- diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
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
    Div[pop_id] <- diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
  }
  
d1 <- c(d1, Div[1])
d2 <- c(d2, Div[2])

}#while

# Create output objects
N_effective <- list()

Offspring_Record1 <-  sapply(1:N, function(x) length(which(Copied[1,] == x)))
Offspring_Record2 <-  sapply(1:N, function(x) length(which(Copied[2,] == x)))

N_effective[[1]] <- N_e(mean(Offspring_Record1), var(Offspring_Record1))
N_effective[[2]] <- N_e(mean(Offspring_Record2), var(Offspring_Record2))


for (t in 1:tmax) {
  
  # Migration
  Migrants1 <- which(rbinom(N,1,m) == 1)
  Migrants2 <- sample(1:N, size = length(Migrants1), replace = FALSE)
  Pop[1, Migrants1] <- Pop[2, Migrants2]
  Pop[2, Migrants2] <- Pop[1, Migrants1]
  
  
  for (pop_id in 1:2) {
    
    #Cultural Transmission


      #Cultural Exchange
      Models_ex <- N + which(rbinom(1:N,1,e) == 1)         #We add N to IDs so we don't confuse them with focal population
      Copied <- sample(c(1:N, Models_ex), N, replace = TRUE)
      
      Copied_pop1 <- Copied[Copied <= N]
      Copied_pop2 <- Copied[Copied > N] - N
      
      Pop[pop_id,] <- c( Pop[pop_id,Copied_pop1], Pop[c(1,2)[-pop_id], Copied_pop2] )
      
      Offspring_Record <-  sapply(c(1:N, Models_ex), function(x) length(which(Copied == x)))
      
    # Compute effective population size
    N_effective[[pop_id]] <-  c(N_effective[[pop_id]], N_e(mean(Offspring_Record), var(Offspring_Record)))
    
    
    # Innovation
    Innovators <- rbinom(N,1,mu)
    Pop[pop_id, Innovators == 1] <- (Counter[pop_id] + 1) : (Counter[pop_id] + length(which(Innovators==1)))
    Counter[pop_id] <- max(Pop[pop_id,])
    
  
  }#pop_id
  
  print(t)
}



par(mfrow=c(1,2))
plot(N_effective[[1]], type = "b")
abline(h = mean(N_effective[[1]]))

plot(N_effective[[2]], type = "b")
abline(h = mean(N_effective[[2]]))










