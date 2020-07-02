
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission

N = 1000          # Census popupation size
tmax = 10000       # Number of timesteps / generations
mu = 10e-4        # Innovation rate
MoT = "OneToMany"
k = 100          # Strength of one-to-many transmission, number of cultural models
f = 2          # Conformity exponent

# Initialize population with cultural traits
Pop <- 1:N

# Create numerator for cultural variants
Counter <- max(N)



Trait_Record     <- matrix(NA, nrow = tmax, ncol = N)
Offspring_Record <- matrix(NA, nrow = tmax, ncol = N)

for (t in 1: tmax){
  print(t)
  
  #Vertical Transmission
  if (MoT == "Vertical"){
  Parents <- sample(1:N, N, replace = TRUE)
  Pop <- Pop[Parents]
  
  for (x in 1:N) {
    Offspring_Record[t,x] <- length(which(Parents == x))
  }
  
  }
  
  #One-to-many Transmission
  if (MoT == "OneToMany"){
  Models <- sample(1:N, size = k, replace = FALSE)
  Copied <- sample(Models, N, TRUE)
  
  Pop <- Pop[Copied]
  
  for (x in 1:N) {
   Offspring_Record[t,x] <- length(which(Copied == x))
  }
  }
  
  #Conformity
  if (MoT == "Conformity"){
  Variants <- unique(Pop)
  
  Freq_Variants <- c()
  
  for (x in Variants) {
    Freq_Variants[which(Variants == x)] <- length(which(Pop[Models] == x))
  }
    
  }

  
  
  # Innovation
  Innovators <- rbinom(N,1,mu)
  Pop[Innovators == 1] <- (Counter + 1) : (Counter + length(which(Innovators==1)))
  Counter <- max(Pop)
  
  
  
  #Record current traits
  Trait_Record[t,] <- Pop

}#t



# Calculate effective population sizes

Offspring_Record <- Offspring_Record[-(1:2000),]
k_bar <- apply(Offspring_Record, 1, mean) 
V_k   <- apply(Offspring_Record, 1, var) 

N_e <- function(k_bar,V_k){
  (N * k_bar - 1 ) / ( (V_k / k_bar) + k_bar  -1 )
}

N_eff <- c()
for (x in 1:80000) {
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

