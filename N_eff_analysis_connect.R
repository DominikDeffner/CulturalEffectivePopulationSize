
# Effective population size analysis script
# Connectedness

# Simulation results

load("~/CulturalEffectivePopulationSize/Neff_connect_0809")
seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-3), m = seq(0,0.9,0.1), e = seq(0,0.9,0.1) )

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
    Upper[i,k] <- quantile(a, 0.95)
    Lower[i,k] <- quantile(a, 0.05)
    
  }
}

Mean <- apply(Mean[, -1], 1, mean)
Upper <- apply(Upper[, -1], 1, mean)
Lower <- apply(Lower[, -1], 1, mean)



graphics.off()
png("Connectedness.png", res = 900, height = 11, width = 18, units = "cm")

par(mfrow=c(1,2),
    oma=c(0,2.2,0,0.5),
    mar=c(4,2,2,0.75))

plot(Mean[which(seq$e == 0)], ylim = c(0, 1100), xaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "Migration")
axis(1,at=seq(1,10,1),labels=seq(0,0.9,0.1))

arrows(seq(1,10,1),Lower[which(seq$e == 0)],seq(1,10,1),Upper[which(seq$e == 0)], code=3, lwd=2,col="darkgray", cex=1.3, length=0, angle = 90)

abline(h = 1000, lty = 2)

mtext(expression(italic(m)), side = 1, line = 2.2, cex = 1.4)



plot(Mean[which(seq$m == 0)], ylim = c(0, 1100), xaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "Cultural Exchange")
axis(1,at=seq(1,10,1),labels=seq(0,0.9,0.1))

arrows(seq(1,10,1),Lower[which(seq$m == 0)],seq(1,10,1),Upper[which(seq$m == 0)], code=3, lwd=2,col="darkgray", cex=1.3, length=0, angle = 90)

abline(h = 1000, lty = 2)
mtext(expression(italic(e)), side = 1, line = 2.2, cex = 1.4)


mtext(expression(paste("Effective population size  ", italic(N[e]))), side = 2, line = 1, cex = 1.2, outer = TRUE)


dev.off()
