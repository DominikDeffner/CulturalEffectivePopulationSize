
# Cultural effective population size model based on Wright-Fisher dynamic

###
##
# Alternative migration mechanism where immigration constantly increases census size in a focal population (section 3.3. and appendix)
##
###


#Load Package to calculate Simpson diversity
library(vegan)

#Function for different Ne formulations (equations 2 and 3 in the main text)

Ne_i <- function(N,k,sigma){
  (N*k-1)  /  (k -1 + (sigma/k))
}

Ne_v <- function(N,k,sigma){
  ((N-1)*k)/(sigma/k)
}


seq<-expand.grid(tmax=150,Nsim = 1000, mu = c(1e-4,1e-3, 1e-2, 1e-1), N_mig = c(5, 10, 50) )



#Define simulation function

sim.funct <- function(tmax, Nsim, mu, N_mig){
  
  Combined_list <- list()
  
  for (sim in 1:Nsim) {
    
    N <-  c(100, 10000)
    
    Pop <- list()
    
    # Maximally diverse
    Pop[[1]] <- sample(1:N[1])
    
    #All same
    Pop[[2]] <- rep(N[1]+1, N[2])
    
    # Create numerator for cultural variants in both populations
    Counter <- max(c(Pop[[1]], Pop[[2]]))
    
    # Calculate Simpson diversity in both populations
    Div <- c()
    for (pop_id in 1:2) {
      Div[pop_id] <- vegan::diversity(sapply(unique(Pop[[pop_id]]), function (x) length(which(Pop[[pop_id]] == x))), index = "simpson")
    }
    
    # Burn-in to reach equilibrium
    while(Div[1] > Div[2]){
      
      Copied <- list()
      for (pop_id in 1:2) {
        #Cultural Transmission
        
        Copied[[pop_id]] <- sample(1:N[pop_id], N[pop_id], replace = TRUE)
        Pop[[pop_id]] <- Pop[[pop_id]][Copied[[pop_id]]]
        
        # Innovation
        Innovators <- rbinom(N[pop_id],1,mu)
        Pop[[pop_id]][Innovators == 1] <- (Counter + 1) : (Counter + length(which(Innovators==1)))
        Counter <- max(c(Pop[[1]], Pop[[2]]))
      }#pop_id
      
      # Calculate Simpson diversity in both populations
      Div <- c()
      for (pop_id in 1:2) {
        Div[pop_id] <- vegan::diversity(sapply(unique(Pop[[pop_id]]), function (x) length(which(Pop[[pop_id]] == x))), index = "simpson")
      }
      
    }#while
    
    
    # Create output objects
    N_effective_v <- list()
    N_effective_i <- list()
    Div_Simpson <- list()
    Div_NoTraits <- list()
    
    # Calculate effective sizes and diversity indices
    
    for (pop_id in 1:2) {
      Offspring_Record <-  sapply(1:N[pop_id], function(x) length(which(Copied[[pop_id]] == x)))
      N_effective_v[[pop_id]] <- Ne_v(N[pop_id],mean(Offspring_Record),var(Offspring_Record))
      N_effective_i[[pop_id]] <- Ne_i(N[pop_id],mean(Offspring_Record),var(Offspring_Record))
      
      Div_Simpson[[pop_id]] <- vegan::diversity(sapply(unique(Pop[[pop_id]]), function (x) length(which(Pop[[pop_id]] == x))), index = "simpson")
      Div_NoTraits[[pop_id]] <- length(unique(Pop[[pop_id]]))
    }
    
    
    for (t in 1:tmax) {
      
      # Migration
      Migrants <- sample(1:length(Pop[[2]]), N_mig, replace = FALSE)
      
      Pop[[1]] <- c( Pop[[1]],  Pop[[2]][Migrants]  )
      Pop[[2]] <- Pop[[2]][-Migrants]
      
      
      for (pop_id in 1:2) {
        
        #Cultural Transmission
        Copied <- sample(1:length(Pop[[pop_id]]), length(Pop[[pop_id]]), replace = TRUE)
        Pop[[pop_id]] <- Pop[[pop_id]][Copied]
        Offspring_Record <-  sapply(1:length(Pop[[pop_id]]), function(x) length(which(Copied == x)))
        
        # Compute effective population sizes
        N_effective_v[[pop_id]] <-  c(N_effective_v[[pop_id]], Ne_v(length(Offspring_Record),mean(Offspring_Record), var(Offspring_Record)))
        N_effective_i[[pop_id]] <-  c(N_effective_i[[pop_id]], Ne_i(length(Offspring_Record),mean(Offspring_Record), var(Offspring_Record)))
        
        #Compute diversity indices
        Div_Simpson[[pop_id]] <- c( Div_Simpson[[pop_id]], vegan::diversity(sapply(unique(Pop[[pop_id]]), function (x) length(which(Pop[[pop_id]] == x))), index = "simpson"))
        Div_NoTraits[[pop_id]] <- c( Div_NoTraits[[pop_id]] ,length(unique(Pop[[pop_id]])))
        
        # Innovation
        Innovators <- rbinom(length(Pop[[pop_id]]),1,mu)
        Pop[[pop_id]][Innovators == 1] <- (Counter + 1) : (Counter + length(which(Innovators==1)))
        Counter <- max(c(Pop[[1]], Pop[[2]]))
        
      }#pop_id
    }
    
    Output_list<-list(N_effective_v = N_effective_v,
                      N_effective_i = N_effective_i,
                      Div_Simpson = Div_Simpson,
                      Div_NoTraits = Div_NoTraits) 
    Combined_list[[sim]]<- Output_list
    
  }#sim
  
  return(Combined_list)  
  
} #end simulation function



# pass to mclapply; it makes sense to select as many cores as there are parameter combinations in case you have access to a computer cluster ("mc.cores" argument)

library(parallel)

result <- mclapply(
  1:nrow(seq) ,
  function(i) sim.funct(seq$tmax[i], seq$Nsim[i], seq$mu[i], seq$N_mig[i]),
  mc.cores=1)


