
load("~/CulturalEffectivePopulationSize/Neff_1000sim_2708")

seq<-expand.grid(N=10000, tmax=300,Nsim = 100, mu = c(1e-2,1e-3,1e-4), k = c(0.25,0.5,0.75,1), theta = c(0.5, 0.9, 1, 1.1, 1.5, 3) )
library(scales)
#color stuff
require(RColorBrewer)#load package


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
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
    }
    Mean[i,k] <- mean(a) 
    Upper[i,k] <- quantile(a,0.95)
    Lower[i,k] <- quantile(a, 0.05)
    
  }
}


graphics.off()
png("Conf.png", res = 600, height = 16, width = 20, units = "cm")

par(mfrow = c(3,3),
    oma=c(2,5.2,2,0),
    mar=c(2.5,2.5,0.5,0.75))
for (mu in c(1e-4,1e-3,1e-2)) {
  for (theta in c(0.5, 1, 1.5)) {
    
    
    i <- which(seq$k==1 & seq$mu==mu & seq$theta==theta)
    
    plot(Mean[i,], type = "l", ylim = c(4000, 10500), col="black", lwd=3, xlim = c(1,100), xlab = "", ylab = "")
    polygon(c(1:300,300:1), c(Upper[i,],rev(Lower[i,])),col=alpha("black",alpha = 0.2), border = NA)
    abline(h = 10000, lty = 2)
    #mtext(seq$mu[i])
    #mtext(seq$theta[i], line = 1.1)
    
  }
}

mtext(expression(paste("Effective population size  ", italic(N[e]))), side = 2,outer=TRUE, line = 3.5, cex = 1.2)
mtext("Generation", side = 1, line = 0.5, outer = TRUE, cex = 1.2)
mtext(expression(paste("   Anti-conformist  (", theta, " = 0.5)", "                       ",
                       "Unbiased  (", theta, " = 1)", "                            ",
                       "Conformist  (", theta, " = 1.5)")), side = 3,outer=TRUE, line = 0.2, cex = 0.9)

mtext(expression(paste("    ",mu, " = ",10^-2, "                                ",
                       mu, " = ",10^-3, "                                ",
                       mu, " = ",10^-4)), side = 2,outer=TRUE, line = 0.2, cex = 0.9)

dev.off()






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