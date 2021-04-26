
# Effective population size analysis script
#Frequency dependence


library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values


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
png("Conf.png", res = 1200, height = 16, width = 18, units = "cm")

par(mfrow = c(3,3),
    oma=c(2,5.2,2,0),
    mar=c(2.5,2.5,0.5,0.75))
for (mu in c(1e-4,1e-3,1e-2)) {
  for (theta in c(0.5,1,1.5)) {
    
    
    i <- which(seq$k==1 & seq$mu==mu & seq$theta==theta)
    
    plot(Mean[i,], type = "n", ylim = c(500, 1100), col="black", lwd=3, xlim = c(1,100), xlab = "", ylab = "")
    polygon(c(1:300,300:1), c(Upper[i,],rev(Lower[i,])),col=alpha(col.pal[1],alpha = 0.5), border = NA)
    lines(Mean[i,], type = "l", ylim = c(500, 1100), col="black", lwd=3)
    
    abline(h = 1000, lty = 2)
    #mtext(seq$mu[i])
    #mtext(seq$theta[i], line = 1.1)
    
  }
}

mtext(expression(paste("Effective population size  ", italic(N[e]))), side = 2,outer=TRUE, line = 3.5, cex = 1.2)
mtext("Generation", side = 1, line = 0.5, outer = TRUE, cex = 1.2)
mtext(expression(paste("   Anti-conformist  (", theta, " = 0.5)", "                 ",
                       "Unbiased  (", theta, " = 1)", "                       ",
                       "Conformist  (", theta, " = 1.5)")), side = 3,outer=TRUE, line = 0.2, cex = 0.9)

mtext(expression(paste("    ",mu, " = ",10^-2, "                                ",
                       mu, " = ",10^-3, "                                ",
                       mu, " = ",10^-4)), side = 2,outer=TRUE, line = 0.2, cex = 0.9)




dev.off()

