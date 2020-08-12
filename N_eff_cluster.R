
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission

N_e <- function(N,k_bar,V_k){
  (N * k_bar - 1 ) / ( (V_k / k_bar) + k_bar  -1 )
}


seq<-expand.grid(N=100, tmax=10,Nsim = 10, mu = 1e-3, k = 1, theta = seq(from = 0.9, to = 1.1, by = 0.0025) )

#Define simulation function

sim.funct <- function(N, tmax, Nsim, mu, k, theta){

Combined_list <- list()

for (sim in 1:Nsim) {
  
  # Initialize population with cultural traits
  Pop <- sample(1:N)
  
  # Create numerator for cultural variants
  Counter <- max(N)
  
  # Output vectors
  
  #Frequency_spectrum <- matrix(NA, ncol = N, nrow = tmax)
  k_bar <- c()
  V_k <- c()
  N_eff <- c()
  
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
    P <- Freq_Variants^theta / sum(Freq_Variants^theta)
    
    #Divide probabilities by frequency and sample from individuals
    P_Ind <- P/Freq_Variants
    Copied <- sample(Models, N, replace = TRUE, sapply(Models, function (x) P_Ind[which(Variants == Pop[x])]))
    Pop <- Pop[Copied]
    
    
    #Record mean and variance of cultural offspring for everyone in the parental generation
    Offspring_Record <-c()
    for (x in 1:N) {
      Offspring_Record[x] <- length(which(Copied == x))
    }#x
    
    k_bar[t] <- mean(Offspring_Record)
    V_k[t]   <- var(Offspring_Record)
    N_eff[t] <-  N_e(N, k_bar[t], V_k[t])
    # Innovation
    Innovators <- rbinom(N,1,mu)
    Pop[Innovators == 1] <- (Counter + 1) : (Counter + length(which(Innovators==1)))
    Counter <- max(Pop)
    
    # Frequency spectrum of traits
    # u <- unique(Pop)
    # f <- sapply(u, function(x) length(which(Pop == x)))
    # 
    # for (i in 1:N) {
    #   Frequency_spectrum[t,i] <- length(which(f == i))
    # }
    
    
  }#t
  
  Output_list<-list(k_bar = k_bar,
                    V_k   = V_k,
                    N_eff = N_eff) 
  
  Combined_list[[sim]]<- Output_list
  
}#sim

return(Combined_list)  

} #end simulation function


# pass to mclapply

library(parallel)

result <- mclapply(
  1:nrow(seq) ,
  function(i) sim.funct(seq$N[i], seq$tmax[i], seq$Nsim[i], seq$mu[i], seq$k[i], seq$theta[i]),
  mc.cores=1)
