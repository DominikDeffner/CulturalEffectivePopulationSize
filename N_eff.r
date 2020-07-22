
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission

N = 1000          # Census popupation size
tmax = 1000       # Number of timesteps / generations
mu = 1e-3        # Innovation rate
k = 1          # Strength of one-to-many transmission, number of cultural models
f = 1.005         # Conformity exponent


# Initialize population with cultural traits
Pop <- sample(1:N)

# Create numerator for cultural variants
Counter <- max(N)

Trait_Record     <- matrix(NA, nrow = tmax, ncol = N)
k_bar <- c()
V_k <- c()

for (t in 1: tmax){
  print(t)
  
  #Cultural Transmission
  
  # First sample set of potential models (fraction k of N) from the population
  Models <- sample(1:N, size = k*N, replace = FALSE)
  
  #Vector with unique variants
  Variants <- unique(Pop[Models])
  
  #Frequency of each variant
  Freq_Variants <- c()
  for (x in Variants) {
    Freq_Variants[which(Variants == x)] <- length(which(Pop[Models] == x))
  }
  
  #Probability individuals choose each variant
  P <- Freq_Variants^f / sum(Freq_Variants^f)
  
  # 2 ways to do conformity
  
  # 1. Sample from variants, then randomly choOse one carrier
  
  #Copied <- c()
  #for (i in 1:N) {
    #if ( length(Variants) == 1) {
   #   Copied[i] <- sample(Models, 1)
  # } else {
     # v <- sample(Variants, size = 1, replace = FALSE, P)
    #  Copied[i] <- sample(Models[which(Pop[Models] == v)], 1)
   # }
  #}

  
  #2. Divide probabilities by frequency and sample from individuals (faster)
   P_Ind <- P/Freq_Variants
   Copied <- sample(Models, N, replace = TRUE, sapply(Models, function (x) P_Ind[which(Variants == Pop[x])]))
  
  Pop <- Pop[Copied]
  

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




par(mfrow = c(2,2))




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
plot(N_eff, type = "l") #, ylim = c(0,(N+0.5*N)))
abline(h=mean(N_eff))
abline(h=mean(N), lty=2)


plot(V_k, type="l")
plot(V_k, N_eff)
abline(lm(N_eff~V_k))


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

