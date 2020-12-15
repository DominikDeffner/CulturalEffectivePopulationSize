
# Effective population size analysis script
# Connectedness

# Simulation results

load("~/CulturalEffectivePopulationSize/Neff_connect_migration_1210")
seq<-expand.grid(N=1000, tmax=300,Nsim = 100, mu = c(1e-4), k = 1, m = seq(0,0.9,0.1), e = seq(0,0.9,0.1) )

Mean_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])


for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    b <- c()

    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective_i"]][[1]][k])
      b <- c(b, 2*result[[i]][[j]][["N_effective_v"]][[1]][k])

    }
    Mean_ei[i,k]  <- mean(a) 
    Upper_ei[i,k] <- quantile(a, 0.95)
    Lower_ei[i,k] <- quantile(a, 0.05)
    
    Mean_ev[i,k] <- mean(b) 
    Upper_ev[i,k] <- quantile(b, 0.95)
    Lower_ev[i,k] <- quantile(b, 0.05)
  }
}

Mean_ei <- apply(Mean_ei[, -1], 1, mean)
Upper_ei <- apply(Upper_ei[, -1], 1, mean)
Lower_ei <- apply(Lower_ei[, -1], 1, mean)

Mean_ev <- apply(Mean_ev[, -1], 1, mean)
Upper_ev <- apply(Upper_ev[, -1], 1, mean)
Lower_ev <- apply(Lower_ev[, -1], 1, mean)




library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values




graphics.off()
png("Connectedness.png", res = 1200, height = 12, width = 20, units = "cm")

par(mfrow=c(1,2),
    oma=c(0,2.2,0,0.5),
    mar=c(3.2,2,2,0.75))



plot((1:10)-0.17,Mean_ev[which(seq$e == 0)], ylim = c(0, 2050),xlim = c(1,10.2), type= "n",xaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "Migration")
axis(1,at=seq(1,10,1),labels=seq(0,0.9,0.1))
abline(h = 1000, lty = 2)

arrows(seq(1,10,1)-0.17,Lower_ev[which(seq$e == 0)],seq(1,10,1)-0.17,Upper_ev[which(seq$e == 0)], code=3, lwd=2,col=col.pal[1], cex=1, length=0, angle = 90)
points((1:10)-0.17, Mean_ev[which(seq$e == 0)],pch=18, cex=1.2, col = col.pal[1])


arrows(seq(1,10,1)+0.17,Lower_ei[which(seq$e == 0)],seq(1,10,1)+0.17,Upper_ei[which(seq$e == 0)], code=3, lwd=2,col=col.pal[3], cex=1, length=0, angle = 90)
points((1:10)+0.17, Mean_ei[which(seq$e == 0)],pch=16, cex=1.2, col = col.pal[3])


mtext(expression(italic(m)), side = 1, line = 2.2, cex = 1.4)

legend("top", c("Variance effective size", "Inbreeding effective size"), col = c(col.pal[1], col.pal[3]), pch = c(18,16), bty = "n", cex = 1.2)



plot((1:10)-0.17,Mean_ev[which(seq$m == 0)], ylim = c(0, 2050),xlim = c(1,10.2), type= "n",xaxt = "n", ylab = "", pch=18, cex=1, xlab = "", main = "Cultural Exchange")
axis(1,at=seq(1,10,1),labels=seq(0,0.9,0.1))
abline(h = 1000, lty = 2)

arrows(seq(1,10,1)-0.17,Lower_ev[which(seq$m == 0)],seq(1,10,1)-0.17,Upper_ev[which(seq$m == 0)], code=3, lwd=2,col=col.pal[1], cex=1, length=0, angle = 90)
points((1:10)-0.17, Mean_ev[which(seq$m == 0)],pch=18, cex=1.2,col = col.pal[1])


arrows(seq(1,10,1)+0.17,Lower_ei[which(seq$m == 0)],seq(1,10,1)+0.17,Upper_ei[which(seq$m == 0)], code=3, lwd=2,col=col.pal[3], cex=1, length=0, angle = 90)
points((1:10)+0.17, Mean_ei[which(seq$m == 0)],pch=16, cex=1.2, col = col.pal[3])

mtext(expression(italic(e)), side = 1, line = 2.2, cex = 1.4)


mtext(expression(paste("Effective population size  ", italic(paste(N[e])))), side = 2, line = 1, cex = 1.2, outer = TRUE)


dev.off()
