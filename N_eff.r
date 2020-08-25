
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission

# We simulate 2 different populations, one starting with all the same variants, one with all unique variants.
# We measure diversity and run simulations until diversities cross indicating they've reached equilibrium, then switch transmission mode

# Function to calculate effective population sizes
library(vegan)

N_e <- function(k_bar,V_k){
  (N * k_bar - 1 ) / ( (V_k / k_bar) + k_bar  -1 )
}

N = 10000          # Census population size
tmax = 300       # Number of timesteps / generations
mu = 1e-4        # Innovation rate
k = 0.5          # Strength of one-to-many transmission, number of cultural models
theta = 1         # Conformity exponent
m = 0.1         #Migration rate between 2 populations

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


for (t in 1:tmax) {
  
   # Migration
  Migrants1 <- rbinom(N,1,m)
  Migrant_Variants <- Pop[1, which(Migrants1 == 1)]
  
  Migrants2 <- sample(1:N, size = length(which(Migrants1==1)), replace = FALSE)

  Pop[1, which(Migrants1 == 1)] <- Pop[2, Migrants2]
  Pop[2, Migrants2] <- Migrant_Variants
  
  for (pop_id in 1:2) {
    
    #Cultural Transmission
    # First sample set of potential models (fraction k of N) from the population
    Models <- sample(1:N, size = k*N, replace = FALSE)
    
    #Vector with unique variants
    Variants <- unique(Pop[pop_id,Models])
    
    #Frequency of each variant
    Freq_Variants <- c()
    for (x in Variants) {
      Freq_Variants[which(Variants == x)] <- length(which(Pop[pop_id,Models] == x))
    }
    
    #Probability individuals choose each variant
    P <- Freq_Variants^theta / sum(Freq_Variants^theta)
    P_Ind <- P/Freq_Variants
    Copied <- sample(Models, N, replace = TRUE, sapply(Models, function (x) P_Ind[which(Variants == Pop[pop_id,x])]))
    
    Pop[pop_id,] <- Pop[pop_id,Copied]
    
    
    # Compute effective population size
    Offspring_Record <-  sapply(1:N, function(x) length(which(Copied == x)))
    N_effective[[pop_id]] <-  c(N_effective[[pop_id]], N_e(mean(Offspring_Record), var(Offspring_Record)))
    
    
    # Innovation
    Innovators <- rbinom(N,1,mu)
    Pop[pop_id, Innovators == 1] <- (Counter[pop_id] + 1) : (Counter[pop_id] + length(which(Innovators==1)))
    Counter[pop_id] <- max(Pop[pop_id,])
    
    
    
    # Frequency spectrum of traits
    
    # Unique Traits
    #u <- unique(Pop)
    
    # Frequency of traits
    #f <- sapply(u, function(x) length(which(Pop == x)))
    #g <- unique(f)
    #z <- sapply(g, function(x) length(which(f == x)))
    
    #for (i in 1:N) {
    #  Frequency_spectra[t,i] <- length(which(f == i))
    #}
    

  }#pop_id
  
  print(t)
}




plot(N_effective[[1]], type = "b")





















par(mfrow = c(2,2))





N_eff <- c()
for (x in 1:length(k_bar)) {
  N_eff[x] <- N_e(k_bar[x], V_k[x])
}
plot(N_eff, type = "l") #, ylim = c(0,(N+0.5*N)))
abline(h=mean(N_eff))
abline(h=mean(N), lty=2)


plot(V_k, type="l")
plot(V_k, N_eff)
abline(lm(N_eff~V_k))

# Trajectories of single traits
Trajectories <- matrix(NA, nrow = tmax, ncol = length(unique(c(Trait_Record))))
for( i in unique(c(Trait_Record))){
  Trajectories[, which(unique(c(Trait_Record))==i)] <- apply(Trait_Record, 1, function(x) length(which(x == i)))
}

Trajectories <- Trajectories / N

a <- Trajectories#[ ,apply(Trajectories, 2, function(x) max(x) > 0.01)]

plot(a[,1], type = "n", ylim = c(0,0.1))
for (j in 1:ncol(a)) {
  lines(a[,j])
}


