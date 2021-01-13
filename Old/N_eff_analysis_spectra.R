
# Effective population size analysis script



#Frequency spectra

seq<-expand.grid(N=10000, tmax=200, Nsim = 100, mu = c(1e-4,1e-3,1e-2,1e-1), k = c(1), theta = seq(from = 0.1, to = 3, by = 0.1) )


N <- unique(seq$N)
theta <- 2

par(mfrow= c(4,6), 
    oma=c(3,2.5,0,0),
    mar=c(1.5,2.4,1.5,0.1))

for (mu in c(1e-4,1e-3,1e-2,1e-1)) {
  
  
  for (t in c(1,5,10, 50, 100, 200)) {
    
    
    #Do some reformatting of the data to make plots look nicer
    
    y <- result[[which(seq$theta == theta & seq$mu == mu)]][[1]]$Frequencies[[t]]
    
    g <- sapply(1:N, function (x) ifelse(y %in% z, return(x), return(0) ) )
    
    z <- result[[which(seq$theta == theta & seq$mu == mu)]][[1]]$Counts[[t]]
    
    z <- sapply(1:N, function (x) ifelse(x %in% g, return(x), return(0) ) )
    
    
    plot(g[order(g)], z[order(g)], pch = 19, xlim = c(1,N),col=alpha("black",alpha = 0.2), log = "xy", ylim = c(1,N), xlab = "Trait Frequency", ylab = "Counts", main = paste("t = ",t))
    
    for (x in 1:  unique(seq$Nsim)) {
      g <- result[[which(seq$theta == theta & seq$mu == mu)]][[x]]$Frequencies[[t]]
      z <- result[[which(seq$theta == theta & seq$mu == mu)]][[x]]$Counts[[t]]
      
      points(g[order(g)], z[order(g)],col=alpha("black",alpha = 0.2), pch=19)
      
    }
    
    mtext(side = 1, "Trait frequency", line = 2, outer = TRUE, cex = 1.2)
    mtext(side = 2, "Counts", line = 0, outer = TRUE, cex = 1.2)
  }
  
}

