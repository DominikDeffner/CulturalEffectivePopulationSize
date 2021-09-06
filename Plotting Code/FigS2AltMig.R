# Simulation results


seq<-expand.grid(tmax=150,Nsim = 1000, mu = c(1e-4,1e-3, 1e-2, 1e-1), N_mig = c(5, 10, 50) )


Mean_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower_ei <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower_ev <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_Simp <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper_Simp <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower_Simp <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_Traits <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper_Traits <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower_Traits <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

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
    Upper_Simp[i,k] <- quantile(c, 0.95)
    Lower_Simp[i,k] <- quantile(c, 0.05)
    
    
    Mean_Traits[i,k]  <- mean(d) 
    Upper_Traits[i,k] <- quantile(d, 0.95)
    Lower_Traits[i,k] <- quantile(d, 0.05)
  }
}

library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values

graphics.off()
png("AltMig.png", res = 600, height = 22, width = 20, units = "cm")

par(mfrow = c(4, 3),
    oma=c(3,3.8,0,0),
    mar=c(2,4,1,1))

for (i in unique(seq$mu)) {
  
Parameters <- which(seq$mu == i & seq$N_mig == 50)

plot(Mean_ei[Parameters,], type = "l", ylim = c(0,8000), col = col.pal[1], lwd = 2, ylab = "Effective population size")

polygon(c(1:150,150:1), c(Upper_ei[Parameters,], rev(Lower_ei[Parameters,])), col=alpha(col.pal[1],alpha = 0.2), border = NA, ylim=c(0,10000))

if (i == 10^-4) mtext(expression(paste(mu, " = ", 10^-4 )), side = 2, line = 6)
if (i == 10^-3) mtext(expression(paste(mu, " = ", 10^-3 )), side = 2, line = 6)
if (i == 10^-2) mtext(expression(paste(mu, " = ", 10^-2 )), side = 2, line = 6)
if (i == 10^-1) mtext(expression(paste(mu, " = ", 10^-1 )), side = 2, line = 6)


plot(Mean_Traits[Parameters,], type = "l", log = "y", ylim = c(1, 10000),lwd = 2, col = col.pal[2], ylab = "Number of unique variants")
polygon(c(1:150,150:1), c(Upper_Traits[Parameters,], rev(Lower_Traits[Parameters,])), col=alpha(col.pal[2],alpha = 0.2), border = NA, ylim=c(0,1000))


plot(Mean_Simp[Parameters,], type = "l", ylim = c(0, 1),lwd = 2, col = col.pal[3], ylab = "Simpson Diversity Index")
polygon(c(1:150,150:1), c(Upper_Simp[Parameters,], rev(Lower_Simp[Parameters,])), col=alpha(col.pal[3],alpha = 0.2), border = NA, ylim=c(0,1))

}
mtext("Generation", side = 1, cex = 1.2, outer = TRUE, line = 2)

dev.off()