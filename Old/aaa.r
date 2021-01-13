
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission

N = 100          # Census popupation size
tmax = 10       # Number of timesteps / generations
mu = 10e-4        # Innovation rate
k = 1          # Strength of one-to-many transmission, number of cultural models
f = 2          # Conformity exponent

# Initialize population with cultural traits
Pop <- 1:N

# Create numerator for cultural variants
Counter <- max(N)

Trait_Record     <- matrix(NA, nrow = tmax, ncol = N)
k_bar <- c()
V_k <- c()

for (t in 1: tmax){
  print(t)
  
  #Cultural Transmission
  Models <- sample(1:N, size = k*N, replace = FALSE)
  Variants <- unique(Pop[Models])
  
  Freq_Variants <- c()
  for (x in Variants) {
    Freq_Variants[which(Variants == x)] <- length(which(Pop[Models] == x))
  }
  
  P <- Freq_Variants^f / sum(Freq_Variants^f)
  
  P_Models <- c()
  
  for (x in Models) {
    P_Models[which(Models == x)] <- P[which(Variants == Pop[x])]
  }

  # Copied <- sample(Models, N, TRUE)
  
  New <- c()
  for (id in 1:N) {
    New[id] <- sample(Variants, 1, prob =  P)
  } 
  
  #Copied <- sample(Variants, N, TRUE, P)
  
  Pop <- New
  
  
  
  
  
  
  #Record mean and variance of cultural offspring for everyone in the parental generation
  Offspring_Record <-c()
  for (x in 1:N) {
    Offspring_Record[x] <- length(which(Copied == x))
  }
  
  k_bar[t] <- mean(Offspring_Record)
  V_k[t]   <- var(Offspring_Record)
  
  # Innovation
  Innovators <- rbinom(N,1,mu)
  Pop[Innovators == 1] <- (Counter + 1) : (Counter + length(which(Innovators==1)))
  Counter <- max(Pop)
  
  #Record current traits
  Trait_Record[t,] <- Pop
  
}#t



# Calculate effective population sizes

k_bar <- k_bar[-(1:(tmax/5))]
V_k   <- V_k[-(1:(tmax/5))]

N_e <- function(k_bar,V_k){
  (N * k_bar - 1 ) / ( (V_k / k_bar) + k_bar  -1 )
}

N_eff <- c()
for (x in 1:length(k_bar)) {
  N_eff[x] <- N_e(k_bar[x], V_k[x])
}
plot(N_eff, type = "l")
abline(h=mean(N_eff))







Trajectories <- matrix(NA, nrow = tmax, ncol = length(unique(c(Trait_Record))))
for( i in unique(c(Trait_Record))){
  Trajectories[, which(unique(c(Trait_Record))==i)] <- apply(Trait_Record, 1, function(x) length(which(x == i)))
}

Trajectories <- Trajectories / N

a <- Trajectories[ ,apply(Trajectories, 2, function(x) max(x) > 0.1)]

plot(a[,1], type = "n", ylim = c(0,1))
for (j in 1:ncol(a)) {
  lines(a[,j])
}

