
# Cultural effective population size model
# We want a simple infinite allele Wright-Fisher-style model with different modes of transmission



library(vegan)
library(igraph)
seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-1,1e-2,1e-3,1e-4),p = seq(0.1,1,0.1),pi = seq(0,1.5,length.out = 10), K = seq(1,10,1),
                 type = c("random","scalefree","smallworld"))

seq <- seq[c(which(seq$type == "random"     & seq$pi == 0  & seq$K ==1),
             which(seq$type == "scalefree"  & seq$p == 0.1 & seq$K ==1),
             which(seq$type == "smallworld" & seq$pi == 0  & seq$p ==0.1)) ,]

#Define simulation function

sim.funct <- function(N, tmax, Nsim, mu, p, pi, K, type){
  
  
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
    
    #Calculate Simpson diversity in both populations
    Div <- c()
    for (pop_id in 1:2) {
      Div[pop_id] <- vegan::diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
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
      for (pop_id in 1:2) {
        Div[pop_id] <- vegan::diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
      }
      
      d1 <- c(d1, Div[1])
      d2 <- c(d2, Div[2])
      
    }#while
    
    
    # Create output objects
    Ne <- function(V_k){
      (N-1) / V_k
    }
    
    N_effective <- list()
    Div_Simpson <- list()
    Div_NoTraits <- list()
    
    
    # Calculate effective sizes and diversity indices
    for (pop_id in 1:2) {
      Offspring_Record <-  sapply(1:N, function(x) length(which(Copied[pop_id,] == x)))
      N_effective[[pop_id]] <- Ne(var(Offspring_Record))
      Div_Simpson[[pop_id]] <- vegan::diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson")
      Div_NoTraits[[pop_id]] <- length(unique(Pop[pop_id,]))
    }
    
    
    if (type == "random"){
      # Simulate random notwork
      g <- erdos.renyi.game(N, p, type = "gnp")
    }
    
    if (type == "scalefree"){
      # Scale free network
      g <- sample_pa(N, power = pi, m = 2, out.dist = NULL, out.seq = NULL, out.pref = FALSE, zero.appeal = 1, directed = FALSE, algorithm ="psumtree", start.graph = NULL)
    }
    
    if (type == "smallworld"){
      # Small world
      g <- watts.strogatz.game(1, N, K, p = 0.01, loops = FALSE, multiple = FALSE)
    }
    
    
    A <- as.matrix(get.adjacency(g, type = "both"))
    
    
    for (t in 1:tmax) {
      
      for (pop_id in 1:2) {
        
        Pop_new <- c()
        Copied <- c()
        
        #Cultural Transmission in network
        for (i in 1:N) {
          if (length(which(A[i,] == 1)) > 1){
            Chosen <- sample(which(A[i,] == 1), 1 )
          } else if (length(which(A[i,] == 1)) == 1){
            Chosen <- which(A[i,] == 1)
          } else {
            Chosen <- i
          }
          Copied <- c(Copied, Chosen)
          Pop_new[i] <- Pop[pop_id, Chosen]  
        }
        
        # Replace population with juveniles
        Pop[pop_id , ] <- Pop_new
        
        # Compute effective population size
        Offspring_Record <-  sapply(1:N, function(x) length(which(Copied == x)))
        N_effective[[pop_id]] <-  c(N_effective[[pop_id]], Ne(var(Offspring_Record)))
        
        
        #Compute diversity indices
        Div_Simpson[[pop_id]] <- c( Div_Simpson[[pop_id]], vegan::diversity(sapply(unique(Pop[pop_id,]), function (x) length(which(Pop[pop_id,] == x))), index = "simpson"))
        Div_NoTraits[[pop_id]] <- c( Div_NoTraits[[pop_id]] ,length(unique(Pop[pop_id,])))
        
        # Innovation
        Innovators <- rbinom(N,1,mu)
        Pop[pop_id, Innovators == 1] <- (Counter[pop_id] + 1) : (Counter[pop_id] + length(which(Innovators==1)))
        Counter[pop_id] <- max(Pop[pop_id,])
        
        
      }#pop_id
      
      print(t)
    }
    
    
    Output_list<-list(N_effective = N_effective,
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
  function(i) sim.funct(seq$N[i], seq$tmax[i], seq$Nsim[i], seq$mu[i], seq$p[i], seq$pi[i],seq$K[i],seq$type[i]),
  mc.cores=60)

save(result, file = "Neff_networks2704newscale")
