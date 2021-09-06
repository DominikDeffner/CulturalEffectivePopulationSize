

# OTM/Diversity Plots

N <- 1000
sigma_c <- function(R) {
  (1000-1)/R 
}


seq<-expand.grid(N=1000, tmax=1000,Nsim = 1000, mu = c(1e-1,1e-2,1e-3,1e-4), k = seq(0.1,1,0.1), theta = 1 )
library(scales)
#color stuff
require(RColorBrewer)#load package


x <- seq(from=0, to=1, by=0.2) # fake data
col.pal <- brewer.pal(length(x), "Dark2") #create a pallette which you loop over for corresponding values



MeanNe <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
UpperNe <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
LowerNe <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

MeanSimp <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
UpperSimp <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
LowerSimp <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

MeanTrait <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
UpperTrait <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
LowerTrait <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    b <- c()
    c <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
      b <- c(b, result[[i]][[j]][["Div_Simpson"]][[1]][k])
      c <- c(c, result[[i]][[j]][["Div_NoTraits"]][[1]][k])
    }
    MeanNe[i,k] <- mean(a) 
    UpperNe[i,k] <- quantile(a,0.95)
    LowerNe[i,k] <- quantile(a, 0.05)
    
    MeanSimp[i,k] <- mean(b) 
    UpperSimp[i,k] <- quantile(b,0.95)
    LowerSimp[i,k] <- quantile(b, 0.05)
    
    MeanTrait[i,k] <- mean(c) 
    UpperTrait[i,k] <- quantile(c,0.95)
    LowerTrait[i,k] <- quantile(c, 0.05)
    
  }
}


MeanNe <-  apply(MeanNe[,100:1000], 1, mean)
LowerNe <- apply(LowerNe[,100:1000], 1, mean)
UpperNe <- apply(UpperNe[,100:1000], 1, mean)


MeanSimp <-  apply(MeanSimp[,100:1000], 1, mean)
LowerSimp <- apply(LowerSimp[,100:1000], 1, mean)
UpperSimp <- apply(UpperSimp[, 100:1000], 1, mean)

MeanTrait <-  apply(MeanTrait[,100:1000], 1, mean)
LowerTrait <- apply(LowerTrait[,100:1000], 1, mean)
UpperTrait <- apply(UpperTrait[,100:1000], 1, mean)

mu <- 1e-4









graphics.off()
png("OTM.png", res = 900, height = 16, width = 18, units = "cm")

par(mfrow=c(2,2),
    oma=c(1.1,0.5,0,0),
    mar=c(2.5,3.5,0.5,0.75))

# Offspring variance 
curve(sigma_c(x), from = 1, to = N, lwd = 3, ylab = "", xlim = c(1,N), log="y",yaxt='n')
mtext(expression(paste("Variance in cultural influence  ",italic(sigma[c]^2))), side = 2, line = 2.3, cex = 1)
legend("topleft", "A", cex=1.1, bty="n")
axis(2,at=c(1,10,100,1000),labels=c(1,10,100,1000))
legend("topright",title = expression(paste("Innovation rate ", mu)), ncol = 2, legend=expression(10^-1,10^-2,10^-3,10^-4), col=c(col.pal[1:4]), pch = c(15,16,17,18), bty = "n", cex = 1.1)

# Effective number
#Analytical predictions
curve( N/ sigma_c(R = x), from = 1, to = N, lwd = 3, ylab = "", xlim = c(1,N), ylim = c(0,1050))
mtext(expression(paste("Effective population size  ", italic(N[e]))), side = 2, line = 2.3, cex = 1)
mtext(expression(paste("Number of role models ", italic(R), " (with ",italic(N), " = 1000)")), side = 1, line = 0, outer = TRUE, cex = 1)
legend("topleft", "B", cex=1.1, bty="n")

for (i in c(400, 800)) {
  segments(x0 = i, y0 = 0, x1 = i, y1 = N/sigma_c(R=i), lty = 2)
  segments(x0 = 0, y0 = i, x1 = i, y1 = N/sigma_c(R=i), lty = 2)
  
}

points(seq(100,1000,100), MeanNe[which(seq$mu == mu)], col = "grey", pch=18, cex=1.5)
arrows(seq(100,1000,100),LowerNe[which(seq$mu == mu)],seq(100,1000,100),UpperNe[which(seq$mu == mu)], code=3, lwd=2,col = "grey", cex=1.3, length=0, angle = 90)




plot(NULL, type = "n", ylim = c(0, 1),xlim = c(50,1000), col="black", lwd=3, xlab = "", ylab = "")
legend("topleft", "C", cex=1.1, bty="n")

for (mu in unique(seq$mu)) {
  i <- which(seq$mu==mu)
  
  lines(MeanNe[i],MeanSimp[i],col=col.pal[which(unique(seq$mu)==mu)], lwd=1.5, type = "b", pch=14+which(unique(seq$mu)==mu))
  
}  
mtext("Simpson Diversity Index", side = 2,outer=FALSE, line = 2.3, cex = 1)


plot(NULL, type = "n", ylim = c(1, 1000),xlim = c(50,1000), log = "y", yaxt = "n",col="black", lwd=3, xlab = "", ylab = "")
axis(2,at=c(1,10,100,1000),labels=c(1,10,100,1000))
legend("topleft", "D", cex=1.1, bty="n")

for (mu in unique(seq$mu)) {
  i <- which(seq$mu==mu)
  lines(MeanNe[i],MeanTrait[i],col=col.pal[which(unique(seq$mu)==mu)],  lwd=1.5, type = "b", pch=14+which(unique(seq$mu)==mu))
}

mtext("Number of unique variants", side = 2,outer=FALSE, line = 2.3, cex = 1)


dev.off()










