
seq<-expand.grid(N=10000, tmax=300,Nsim = 100, mu = c(1e-4), k = 1, theta = seq(0.1,3,0.05) )
library(scales)
#color stuff
require(RColorBrewer)#load package



# N_eff curves over time

x <- seq(from=0, to=1, by=0.2) # fake data
col.pal <- brewer.pal(length(x), "Dark2") #create a pallette which you loop over for corresponding values


Mean <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[2]][k])
    }
    Mean[i,k] <- mean(a) 
    Upper[i,k] <- quantile(a, 1)
    Lower[i,k] <- quantile(a, 0)
    
  }
}
  
par(mfrow= c(1,5))

for (i in which(seq$theta %in% c(0.4, 0.9, 1, 1.1, 3) )) {
  
  if (seq$theta[i] %in% c(0.9,1,1.1)){
    x <- c(9500, 10500)
  } else {
    x <- c(3000, 10500)
  }
  
  plot(Mean[i,], type = "l", ylim = x, xlim=c(0,100), col=col.pal[1], lwd=3)
  polygon(c(1:300,300:1), c(Upper[i,],rev(Lower[i,])),col=alpha(col.pal[1],alpha = 0.2), border = NA)
  abline(h = 10000)
  mtext(seq$mu[i])
  mtext(seq$theta[i], line = 1.1)
  
}



#Frequency spectra

seq<-expand.grid(N=10000, tmax=300,Nsim = 100, mu = c(1e-4), k = 1, theta = seq(0.1,3,0.05) )


pop <- 2
mu <- 1e-4

par(mfrow= c(5,10), 
    oma=c(3,2.5,0,0),
    mar=c(1.5,2.4,1.5,0.1))

for (theta in c(0.1, 0.9, 1, 1.1, 3)) {
  
  
  for (t in c(1,2,3,4,5,10,20, 50, 100,300)) {
    
    
    #Do some reformatting of the data to make plots look nicer
    
    #Frequencies
    g   <- c(result[[which(seq$theta == theta)]][[1]]$Frequencies[[pop]][[t]])
    
    # Counts of frequencies
    z <- c(result[[which(seq$theta == theta)]][[1]]$Counts[[pop]][[t]])
    
    
    #g <- sapply(1:N, function (x) ifelse(x %in% g, return(x), return(0) ) )
    
    
   # z <- sapply(1:N, function (x) ifelse(x %in% g, return(x), return(0) ) )
    
    
    plot(g[order(g)], z[order(g)], pch = 19, xlim = c(1,N),col=alpha(col.pal[1],alpha = 0.2), log = "xy", ylim = c(1,N), xlab = "Trait Frequency", ylab = "Counts", main = paste("t = ",t))
    
    for (x in 1:  unique(seq$Nsim)) {
      g <- result[[which(seq$theta == theta)]][[x]]$Frequencies[[pop]][[t]]
      z <- result[[which(seq$theta == theta)]][[x]]$Counts[[pop]][[t]]
      
      points(g[order(g)], z[order(g)],col=alpha(col.pal[1],alpha = 0.2), pch=19)
      
    }
    
    mtext(side = 1, "Trait frequency", line = 2, outer = TRUE, cex = 1.2)
    mtext(side = 2, "Counts", line = 0, outer = TRUE, cex = 1.2)
  }
  
}