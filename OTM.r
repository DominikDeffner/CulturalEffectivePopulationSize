
 # This is Laurel's OTM wizardry

N = 1000

sigma_c <- function(N = n, R) {
  (N * (1 - 1/R) + ( ( R * (N - R)) / N ) * ((N/R) - 0)^2 ) / N
}

# Simulation results

seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-3), k = seq(0.1,1,0.1), theta = 1 )
Mean <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective_s"]][[1]][k])
    }
    Mean[i,k] <- mean(a) 
    Upper[i,k] <- quantile(a, 0.95)
    Lower[i,k] <- quantile(a, 0.05)
    
  }
}

Mean <- apply(Mean[, -1], 1, mean)
Upper <- apply(Upper[, -1], 1, mean)
Lower <- apply(Lower[, -1], 1, mean)



graphics.off()
png("OTM.png", res = 900, height = 9, width = 18, units = "cm")

par(mfrow=c(1,2),
    oma=c(1.1,0,0,0),
    mar=c(2.5,3.5,0.5,0.75))

# Offspring variance 
curve(sigma_c(R = x), from = 1, to = N, lwd = 3, ylab = "", xlim = c(1,n), log="y",yaxt='n')
mtext(expression(paste("Variance in learner number  ",italic(sigma[c]^2))), side = 2, line = 2.2, cex = 0.9)

axis(2,at=c(1,10,100,1000),labels=c(1,10,100,1000))

# Effective number
#Analytical predictions
curve(N / sigma_c(R = x), from = 1, to = N, lwd = 3, ylab = "", xlim = c(1,n), ylim = c(0,1050))
mtext(expression(paste("Effective population size  ", italic(N[e]))), side = 2, line = 2.2, cex = 0.9)
mtext(expression(paste("Number of teaching individuals ", italic(R), " (with ",italic(N), " = 1000)")), side = 1, line = 0, outer = TRUE, cex = 0.9)


for (i in c(400, 800)) {
  segments(x0 = i, y0 = 0, x1 = i, y1 = N/sigma_c(R=i), lty = 2)
  segments(x0 = 0, y0 = i, x1 = i, y1 = N/sigma_c(R=i), lty = 2)
  
}

points(seq(100,1000,100), Mean, col = "darkgray", pch=18, cex=1.3)
arrows(seq(100,1000,100),Lower,seq(100,1000,100),Upper, code=3, lwd=2,col="darkgray", cex=1.3, length=0, angle = 90)


#Simulation results

dev.off()
