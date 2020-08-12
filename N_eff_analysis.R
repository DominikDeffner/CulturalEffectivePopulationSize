
seq<-expand.grid(N=10000, tmax=200, Nsim = 100, mu = c(1e-4,1e-3,1e-2,1e-1), k = c(1), theta = seq(from = 0.1, to = 3, by = 0.1) )
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
      a <- c(a, result[[i]][[j]][["N_eff"]][k])
    }
    Mean[i,k] <- mean(a) 
    Upper[i,k] <- quantile(a, 1)
    Lower[i,k] <- quantile(a, 0)
    
  }
}

par(mfrow = c(1,4))
for (i in 1 : nrow(seq)) {
  
  plot(Mean[i,], type = "l", ylim = c(3000, 10500), col=col.pal[1], lwd=3)
  polygon(c(1:200,200:1), c(Upper[i,],rev(Lower[i,])),col=alpha(col.pal[1],alpha = 0.2), border = NA)
  abline(h = 10000)
  mtext(seq$mu[i])
  mtext(seq$theta[i], line = 1.1)
  
}



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
    
    #Frequencies
    g   <- result[[which(seq$theta == theta & seq$mu == mu)]][[1]]$Frequencies[[t]]
    
    # Counts of frequencies
    z <- result[[which(seq$theta == theta & seq$mu == mu)]][[1]]$Counts[[t]]
    
    
    g <- sapply(1:N, function (x) ifelse(x %in% g, return(x), return(0) ) )
    
    
    z <- sapply(1:N, function (x) ifelse(x %in% g, return(x), return(0) ) )
    
    
    plot(g[order(g)], z[order(g)], type="l", pch = 19, xlim = c(1,N),col=alpha(col.pal[1],alpha = 0.2), log = "xy", ylim = c(1,N), xlab = "Trait Frequency", ylab = "Counts", main = paste("t = ",t))
    
    for (x in 1:  unique(seq$Nsim)) {
      g <- result[[which(seq$theta == theta & seq$mu == mu)]][[x]]$Frequencies[[t]]
      z <- result[[which(seq$theta == theta & seq$mu == mu)]][[x]]$Counts[[t]]
      
      lines(g[order(g)], z[order(g)],col=alpha(col.pal[1],alpha = 0.2), pch=19)
      
    }
    
    mtext(side = 1, "Trait frequency", line = 2, outer = TRUE, cex = 1.2)
    mtext(side = 2, "Counts", line = 0, outer = TRUE, cex = 1.2)
  }
  
}