
# Effective population size analysis script
# Connectedness

#color stuff
require(RColorBrewer)#load package

x <- seq(from=0, to=1, by=0.2) # fake data
col.pal <- brewer.pal(length(x), "Dark2") #create a pallette which you loop over for corresponding values

# Simulation results

#seq<-expand.grid(N=1000, tmax=500,Nsim = 10, mu = c(1e-1,1e-2,1e-3,1e-4), m = seq(0,0.9,0.1), e = seq(0,0.9,0.1) )


Mean_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_Simp <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Mean_Traits <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    b <- c()
    c <- c()
    d <- c()
    
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective_i"]][[1]][k])
      b <- c(b, result[[i]][[j]][["N_effective_v"]][[1]][k])
      c <- c(c, result[[i]][[j]][["Div_Simpson"]][[1]][k])
      d <- c(d, result[[i]][[j]][["Div_NoTraits"]][[1]][k])
      
    }
    Mean_ei[i,k]  <- mean(a) 
    Upper_ei[i,k] <- quantile(a, 0.95)
    Lower_ei[i,k] <- quantile(a, 0.05)
    
    Mean_ev[i,k] <- mean(b) 
    Upper_ev[i,k] <- quantile(b, 0.95)
    Lower_ev[i,k] <- quantile(b, 0.05)
    
    Mean_Simp[i,k]  <- mean(c) 
    Mean_Traits[i,k]  <- mean(d) 
    
  }
}

Mean_ei <- apply(Mean_ei[, -1], 1, mean)
Upper_ei <- apply(Upper_ei[, -1], 1, mean)
Lower_ei <- apply(Lower_ei[, -1], 1, mean)

Mean_ev <- apply(Mean_ev[, -1], 1, mean)
Upper_ev <- apply(Upper_ev[, -1], 1, mean)
Lower_ev <- apply(Lower_ev[, -1], 1, mean)

Mean_Simp <- apply(Mean_Simp[, -1], 1, mean)
Mean_Traits <- apply(Mean_Traits[, -1], 1, mean)








graphics.off()
png("Connectedness.png", res = 900, height = 19, width = 15, units = "cm")

par(mfrow=c(3,2),
    oma=c(3,1,0,0.5),
    mar=c(1.5,3.5,2,0))

#Effective Size
mu <- 1e-4
plot((1:10)-0.17,Mean_ev[which(seq$e == 0 & seq$mu == mu)], ylim = c(0, 2050),xlim = c(1,10.2), type= "n",xaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "")
axis(1,at=seq(1,10,3),labels=seq(0,0.9,0.3))
abline(h = 1000, lty = 2)
mtext("Migration", side = 3, line = 0.5, cex = 1.2)

arrows(seq(1,10,1)-0.17,Lower_ev[which(seq$e == 0 & seq$mu == mu)],seq(1,10,1)-0.17,Upper_ev[which(seq$e == 0 & seq$mu == mu)], code=3, lwd=2,col=col.pal[5], cex=1, length=0, angle = 90)
points((1:10)-0.17, Mean_ev[which(seq$e == 0 & seq$mu == mu)],pch=0, cex=1.2, col = col.pal[5])


arrows(seq(1,10,1)+0.17,Lower_ei[which(seq$e == 0 & seq$mu == mu)],seq(1,10,1)+0.17,Upper_ei[which(seq$e == 0 & seq$mu == mu)], code=3, lwd=2,col=col.pal[6], cex=1, length=0, angle = 90)
points((1:10)+0.17, Mean_ei[which(seq$e == 0 & seq$mu == mu)],pch=1, cex=1.2, col = col.pal[6])


mtext(expression(paste("Effective population size  ", italic(paste(N[e])))), side = 2, line = 3, cex = 1)



plot((1:10)-0.17,Mean_ev[which(seq$m == 0 & seq$mu == mu)], ylim = c(0, 2050),xlim = c(1,10.2), type= "n",xaxt = "n", ylab = "", pch=18, cex=1, xlab = "")
axis(1,at=seq(1,10,3),labels=seq(0,0.9,0.3))
abline(h = 1000, lty = 2)
mtext("Cultural Exchange", side = 3, line = 0.5, cex = 1.2)

arrows(seq(1,10,1)-0.17,Lower_ev[which(seq$m == 0 & seq$mu == mu)],seq(1,10,1)-0.17,Upper_ev[which(seq$m == 0 & seq$mu == mu)], code=3, lwd=2,col=col.pal[5], cex=1, length=0, angle = 90)
points((1:10)-0.17, Mean_ev[which(seq$m == 0 & seq$mu == mu)],pch=0, cex=1.2,col = col.pal[5])

arrows(seq(1,10,1)+0.17,Lower_ei[which(seq$m == 0 & seq$mu == mu)],seq(1,10,1)+0.17,Upper_ei[which(seq$m == 0 & seq$mu == mu)], code=3, lwd=2,col=col.pal[6], cex=1, length=0, angle = 90)
points((1:10)+0.17, Mean_ei[which(seq$m == 0 & seq$mu == mu)],pch=1, cex=1.2, col = col.pal[6])
legend("bottomright", c("Variance effective size", "Inbreeding effective size"), col = c(col.pal[5], col.pal[6]), pch = c(0,1), bty = "n", cex = 1.3)


###
###
#Simpson
###
###

plot(Mean_Simp[which(seq$e == 0& seq$mu == mu)], ylim = c(0, 1),type = "n", xaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "")
axis(1,at=seq(1,10,3),labels=seq(0,0.9,0.3))

for (mu in  unique(seq$mu)) {
  lines(Mean_Simp[which(seq$e == 0& seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
}

mtext("Simpson Diversity Index", side = 2, line = 3, cex = 1)


plot(Mean_Simp[which(seq$m == 0& seq$mu == mu)], ylim = c(0, 1),type = "n", xaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "")
axis(1,at=seq(1,10,3),labels=seq(0,0.9,0.3))

for (mu in unique(seq$mu)) {
  lines(Mean_Simp[which(seq$m == 0& seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2,col=col.pal[which(unique(seq$mu)==mu)])
}
abline(h = 1000, lty = 2)
legend("bottomright",title = expression(paste("Innovation rate ", mu)), ncol = 2, legend=expression(10^-1,10^-2,10^-3,10^-4), col=c(col.pal[1:4]), pch = c(15,16,17,18), bty = "n", cex = 1.3)

#
#
#
#Traits
#
#
#

plot(Mean_Traits[which(seq$e == 0& seq$mu == mu)],log="y", ylim = c(1, 1000),type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "")
axis(1,at=seq(1,10,3),labels=seq(0,0.9,0.3))
axis(2,at=c(1,10,100,1000),labels=c(1,10,100,1000))

for (mu in  unique(seq$mu)) {
  lines(Mean_Traits[which(seq$e == 0& seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
}

mtext(expression(paste("Migration rate  ", italic(paste(m)))), side = 1, line = 3.5, cex = 1)
mtext("Number of unique variants", side = 2, line = 3, cex = 1)


plot(Mean_Traits[which(seq$m == 0& seq$mu == mu)],log="y", ylim = c(1, 1000),type = "n", xaxt = "n",yaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "")
axis(1,at=seq(1,10,3),labels=seq(0,0.9,0.3))
axis(2,at=c(1,10,100,1000),labels=c(1,10,100,1000))

for (mu in unique(seq$mu)) {
  lines(Mean_Traits[which(seq$m == 0& seq$mu == mu)],type = "b",pch=14+which(unique(seq$mu)==mu), cex=1.2,col=col.pal[which(unique(seq$mu)==mu)])
}
mtext(expression(paste("Exchange rate  ", italic(paste(e)))), side = 1, line = 3.5, cex = 1)

dev.off()

