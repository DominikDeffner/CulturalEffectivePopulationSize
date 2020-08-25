
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission



library(vegan)

seq<-expand.grid(N=10000, tmax=300,Nsim = 100, mu = c(1e-4), k = 1, theta = seq(0.1,3,0.05) )

#Define simulation function

sim.funct <- function(N, tmax, Nsim, mu, k, theta){
  
  
  
  Combined_list <- list()
  
  for (sim in 1:Nsim) {
    
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
    
    
    N_effective <- list()
    Frequencies <- list()
    Frequencies[[1]] <- list()
    Frequencies[[2]] <- list()
    
    Counts      <- list()
    Counts[[1]]      <- list()
    Counts[[2]]      <- list()
    
    Offspring_Record1 <-  sapply(1:N, function(x) length(which(Copied[1,] == x)))
    Offspring_Record2 <-  sapply(1:N, function(x) length(which(Copied[2,] == x)))
    
    N_effective[[1]] <- (N * mean(Offspring_Record1) - 1 ) / ( (var(Offspring_Record1)) + mean(Offspring_Record1)  -1 )
    N_effective[[2]] <- (N * mean(Offspring_Record2) - 1 ) / ( (var(Offspring_Record2)) + mean(Offspring_Record2)  -1 )
    
    
    for (t in 1:tmax) {
      
      for (pop_id in 1:2) {
        
        # Frequency spectrum of traits
        
        # Unique Traits
        u <- unique(Pop[pop_id,])
        
        # Frequency of traits
        f <- sapply(u, function(x) length(which(Pop[pop_id,] == x)))
        g <- unique(f)
        z <- sapply(g, function(x) length(which(f == x)))
        
        Frequencies[[pop_id]][[t]] <- g
        Counts[[pop_id]][[t]]      <- z
        
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
        Ne <- (N * mean(Offspring_Record) - 1 ) / ( (var(Offspring_Record)) + mean(Offspring_Record)  -1 )
        N_effective[[pop_id]] <-  c(N_effective[[pop_id]],Ne)    
        
        
        # Innovation
        Innovators <- rbinom(N,1,mu)
        Pop[pop_id, Innovators == 1] <- (Counter[pop_id] + 1) : (Counter[pop_id] + length(which(Innovators==1)))
        Counter[pop_id] <- max(Pop[pop_id,])
        
        
        
        
      }#pop_id
      
    }
    
    
    Output_list<-list(N_effective = N_effective, 
                      Frequencies = Frequencies, 
                      Counts = Counts) 
    
    Combined_list[[sim]]<- Output_list
    
  }#sim
  
  return(Combined_list)  
  
} #end simulation function



# pass to mclapply

library(parallel)

result <- mclapply(
  1:nrow(seq) ,
  function(i) sim.funct(seq$N[i], seq$tmax[i], seq$Nsim[i], seq$mu[i], seq$k[i], seq$theta[i]),
  mc.cores=59)

save(result, file = "Neff_1908_2popspectra")