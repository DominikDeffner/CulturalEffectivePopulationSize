
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission



library(vegan)

Ne_simpl <- function(N,sigma){
  (N-1)/sigma
}

#For OTM
seq<-expand.grid(N=1000, tmax=300,Nsim = 100, mu = c(1e-2,1e-3,1e-4), k = seq(0.1,1,0.1), theta = 1 )

#For Frequency dependence
#seq<-expand.grid(N=1000, tmax=300,Nsim = 100, mu = c(1e-2,1e-3,1e-4), k = 1, theta = c(0.5,1,1.5) )


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
    Div_Simpson <- list()
    Div_Shannon <- list()
    Div_NoTraits <- list()
    
    # Calculate effective population size
    Offspring_Record1 <-  sapply(1:N, function(x) length(which(Copied[1,] == x)))
    Offspring_Record2 <-  sapply(1:N, function(x) length(which(Copied[2,] == x)))
    
    N_effective[[1]] <- Ne_simpl(N,var(Offspring_Record1))
    N_effective[[2]] <- Ne_simpl(N,var(Offspring_Record2))
    
    # Calculate diversity indices

    for (pop_id in 1:2) {
      Div_Simpson[[pop_id]] <- diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
      Div_Shannon[[pop_id]] <- diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "shannon")
      
      Div_NoTraits[[pop_id]] <- length(unique(Pop[pop_id,]))
    }
    
    
    
    for (t in 1:tmax) {
      
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
        Ne <- Ne_simpl(N, var(Offspring_Record))
        N_effective[[pop_id]] <-  c(N_effective[[pop_id]],Ne)    
        
        #Compute diversity indices
        Div_Simpson[[pop_id]] <- c( Div_Simpson[[pop_id]], diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson"))
        Div_Shannon[[pop_id]] <- c( Div_Shannon[[pop_id]], diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "shannon"))
        
        Div_NoTraits[[pop_id]] <- c( Div_NoTraits[[pop_id]] ,length(unique(Pop[pop_id,])))
        
        
        # Innovation
        Innovators <- rbinom(N,1,mu)
        Pop[pop_id, Innovators == 1] <- (Counter[pop_id] + 1) : (Counter[pop_id] + length(which(Innovators==1)))
        Counter[pop_id] <- max(Pop[pop_id,])
        
        
        
        
      }#pop_id
      
    }
    
    
    Output_list<-list(N_effective = N_effective,
                      Div_Shannon = Div_Shannon,
                      Div_Simpson = Div_Simpson,
                      Div_NoTraits = Div_NoTraits) 
    
    Combined_list[[sim]]<- Output_list
    
  }#sim
  
  return(Combined_list)  
  
} #end simulation function



# pass to mclapply

library(parallel)

result <- mclapply(
  1:nrow(seq) ,
  function(i) sim.funct(seq$N[i], seq$tmax[i], seq$Nsim[i], seq$mu[i], seq$k[i], seq$theta[i]),
  mc.cores=30)

save(result, file = "Neff_DiversityOTM1112")