# Conformity/Diversity Plots

seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-2,1e-3,1e-4), k = 1, theta = c(0.5,1,1.5) )
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
      a <- c(a, result2[[i]][[j]][["Div_Simpson"]][[1]][k])
    }
    Mean[i,k] <- mean(a) 
    Upper[i,k] <- quantile(a,0.95)
    Lower[i,k] <- quantile(a, 0.05)
    
  }
}


graphics.off()
png("ConfSimp.png", res = 600, height = 16, width = 20, units = "cm")

par(mfrow = c(3,3),
    oma=c(2,5.2,2,0),
    mar=c(2.5,2.5,0.5,0.75))
for (mu in c(1e-4,1e-3,1e-2)) {
  for (theta in c(0.5, 1, 1.5)) {
    
    
    i <- which(seq$k==1 & seq$mu==mu & seq$theta==theta)
    
    plot(Mean[i,], type = "l", ylim = c(0, 1), col="black", lwd=3, xlim = c(1,300), xlab = "", ylab = "")
    polygon(c(1:300,300:1), c(Upper[i,],rev(Lower[i,])),col=alpha("black",alpha = 0.2), border = NA)
    abline(h = 1000, lty = 2)
    #mtext(seq$mu[i])
    #mtext(seq$theta[i], line = 1.1)
    
  }
}

mtext("Simpson Diversity Index", side = 2,outer=TRUE, line = 3.5, cex = 1.2)
mtext("Generation", side = 1, line = 0.5, outer = TRUE, cex = 1.2)
mtext(expression(paste("   Anti-conformist  (", theta, " = 0.5)", "                       ",
                       "Unbiased  (", theta, " = 1)", "                            ",
                       "Conformist  (", theta, " = 1.5)")), side = 3,outer=TRUE, line = 0.2, cex = 0.9)

mtext(expression(paste("    ",mu, " = ",10^-2, "                                ",
                       mu, " = ",10^-3, "                                ",
                       mu, " = ",10^-4)), side = 2,outer=TRUE, line = 0.2, cex = 0.9)

dev.off()



Mean <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result2[[i]][[j]][["Div_NoTraits"]][[1]][k])
    }
    Mean[i,k] <- mean(a) 
    Upper[i,k] <- quantile(a,0.95)
    Lower[i,k] <- quantile(a, 0.05)
    
  }
}


graphics.off()
png("ConfTraits.png", res = 600, height = 16, width = 20, units = "cm")

par(mfrow = c(3,3),
    oma=c(2,5.2,2,0),
    mar=c(2.5,2.5,0.5,0.75))
for (mu in c(1e-4,1e-3,1e-2)) {
  for (theta in c(0.5, 1, 1.5)) {
    
    
    i <- which(seq$k==1 & seq$mu==mu & seq$theta==theta)
    
    plot(Mean[i,], type = "l", ylim = c(1, 270), col="black", lwd=3, xlim = c(1,300), xlab = "", ylab = "")
    polygon(c(1:300,300:1), c(Upper[i,],rev(Lower[i,])),col=alpha("black",alpha = 0.2), border = NA)
    abline(h = 1000, lty = 2)
    #mtext(seq$mu[i])
    #mtext(seq$theta[i], line = 1.1)
    
  }
}

mtext("Number of Traits", side = 2,outer=TRUE, line = 3.5, cex = 1.2)
mtext("Generation", side = 1, line = 0.5, outer = TRUE, cex = 1.2)
mtext(expression(paste("   Anti-conformist  (", theta, " = 0.5)", "                       ",
                       "Unbiased  (", theta, " = 1)", "                            ",
                       "Conformist  (", theta, " = 1.5)")), side = 3,outer=TRUE, line = 0.2, cex = 0.9)

mtext(expression(paste("    ",mu, " = ",10^-2, "                                ",
                       mu, " = ",10^-3, "                                ",
                       mu, " = ",10^-4)), side = 2,outer=TRUE, line = 0.2, cex = 0.9)

dev.off()











###
##
#
#OTM
#
##
###


# OTM/Diversity Plots


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


graphics.off()
png("OTMDiv.png", res = 600, height = 9, width = 18, units = "cm")

par(mfrow = c(1,2),
    oma=c(1,0,0,0),
    mar=c(2.5,3.5,0.5,0.75))

plot(NULL, type = "n", ylim = c(0, 1),xlim = c(50,1000), col="black", lwd=3, xlab = "", ylab = "")

for (mu in unique(seq$mu)) {
  i <- which(seq$mu==mu)
  
  lines(MeanNe[i],MeanSimp[i],col=col.pal[which(unique(seq$mu)==mu)], lwd=1.5, type = "b", pch=16)
  #arrows(MeanNe[i],LowerSimp[i],MeanNe[i],UpperSimp[i], code=3, lwd=2,col="darkgray", cex=1.3, length=0, angle = 90)
  
}  
mtext("Simpson Diversity Index", side = 2,outer=FALSE, line = 2.3, cex = 1.2)


plot(NULL, type = "n", ylim = c(0, 330),xlim = c(50,1000), col="black", lwd=3, xlab = "", ylab = "")

for (mu in unique(seq$mu)) {
  i <- which(seq$mu==mu)
  
  lines(MeanNe[i],MeanTrait[i],col=col.pal[which(unique(seq$mu)==mu)],  lwd=1.5, type = "b", pch=16)
  arrows(MeanNe[i],LowerTrait[i],MeanNe[i],UpperTrait[i], code=3, lwd=2,col=col.pal[which(unique(seq$mu)==mu)], cex=1.3, length=0, angle = 90)
  
}

mtext(expression(paste("Effective population size  ", italic(N[e]))), side = 1,outer=TRUE, line = 0, cex = 1.2)
mtext("Number of Traits", side = 2,outer=FALSE, line = 2.3, cex = 1.2)

legend("topleft",title = expression(paste("Innovation rate ", mu)), legend=expression(10^-1,10^-2,10^-3,10^-4), col=c(col.pal[1:4]), lty = 1, lwd=2, bty = "n", cex = 0.8)

dev.off()










