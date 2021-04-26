
# Cultural effective population size model
# Connected model; every new innovation is unique (has not ocurred in any population)

library(vegan)

Ne_i <- function(N,k,sigma){
  (N*k-1)  /  (k -1 + (sigma/k))
}

Ne_v <- function(N,k,sigma){
  ((N-1)*k)/(sigma/k)
}

seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-4,1e-3,1e-2,1e-1), m = seq(0,0.9,0.1), e = seq(0,0.9,0.1) )
seq <- seq[-which(seq$m > 0 & seq$e !=0), ]

#Define simulation function

sim.funct <- function(N, tmax, Nsim, mu, m, e){
  
  Combined_list <- list()
  
  for (sim in 1:Nsim) {
    
    # Initialize population with cultural traits
    Pop <- matrix(NA, 2, N)
    
    # Maximally diverse
    Pop[1,] <- sample(1:N)
    
    #All same
    Pop[2,] <- rep(N+1, N) #Make variants in Population 2 unique
    
    # Create numerator for cultural variants in both populations
    Counter <- max(Pop)
    
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
        Pop[pop_id, Innovators == 1] <- (Counter + 1) : (Counter + length(which(Innovators==1)))
        Counter <- max(Pop)
      }#pop_id
      
      # Calculate Simpson diversity in both populations
      Div <- c()
      for (pop_id in 1:2) {
        Div[pop_id] <- diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
      }
    }#while
    
    
    # Create output objects
    N_effective_v <- list()
    N_effective_i <- list()
    Div_Simpson <- list()
    Div_NoTraits <- list()
    
    # Calculate effective sizes and diversity indices
    
    for (pop_id in 1:2) {
      Offspring_Record <-  sapply(1:N, function(x) length(which(Copied[pop_id,] == x)))
      N_effective_v[[pop_id]] <- Ne_v(N,mean(Offspring_Record),var(Offspring_Record))
      N_effective_i[[pop_id]] <- Ne_i(N,mean(Offspring_Record),var(Offspring_Record))
      
      Div_Simpson[[pop_id]] <- diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
      Div_NoTraits[[pop_id]] <- length(unique(Pop[pop_id,]))
    }
    
    
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
        
        Copied_pop_own <- Copied[Copied <= N]
        Copied_pop_other <- Copied[Copied > N] - N
        
        Pop[pop_id,] <- c( Pop[pop_id,Copied_pop_own], Pop[c(1,2)[-pop_id], Copied_pop_other] )
        
        Offspring_Record <-  sapply(c(1:N, Models_ex), function(x) length(which(Copied == x)))
        
        # Compute effective population sizes
        N_effective_v[[pop_id]] <-  c(N_effective_v[[pop_id]], Ne_v(length(Offspring_Record),mean(Offspring_Record), var(Offspring_Record)))
        N_effective_i[[pop_id]] <-  c(N_effective_i[[pop_id]], Ne_i(length(Offspring_Record),mean(Offspring_Record), var(Offspring_Record)))
        
        #Compute diversity indices
        Div_Simpson[[pop_id]] <- c( Div_Simpson[[pop_id]], diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson"))
        Div_NoTraits[[pop_id]] <- c( Div_NoTraits[[pop_id]] ,length(unique(Pop[pop_id,])))
        
        # Innovation
        Innovators <- rbinom(N,1,mu)
        Pop[pop_id, Innovators == 1] <- (Counter + 1) : (Counter + length(which(Innovators==1)))
        Counter <- max(Pop)
        
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



# pass to mclapply

library(parallel)

result <- mclapply(
  1:nrow(seq) ,
  function(i) sim.funct(seq$N[i], seq$tmax[i], seq$Nsim[i], seq$mu[i], seq$m[i], seq$e[i]),
  mc.cores=76)

save(result, file = "Neff_Conn_unique200421")