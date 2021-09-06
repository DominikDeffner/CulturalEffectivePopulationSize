####
###
##
# Plotting code for Fig. 5
##
###
####

# Requires simulation results from "SimNetworks.R" 

library(igraph)

seq<-expand.grid(N=1000, tmax=300,Nsim = 100, mu = c(1e-1,1e-2,1e-3,1e-4),p = seq(0.1,1,0.1),pi = seq(0,1.5,length.out = 10), K = seq(1,10,1),
                 type = c("random","scalefree","smallworld"))

seq <- seq[c(which(seq$type == "random"     & seq$pi == 0  & seq$K ==1),
             which(seq$type == "scalefree"  & seq$p == 0.1 & seq$K ==1),
             which(seq$type == "smallworld" & seq$pi == 0  & seq$p ==0.1)) ,]

Mean <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Upper <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Lower <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_Traits <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
Mean_Simp <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])


for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    b <- c()
    d <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
      b <- c(b, result[[i]][[j]][["Div_NoTraits"]][[1]][k])
      d <- c(d, result[[i]][[j]][["Div_Simpson"]][[1]][k])
      
    }
    Mean[i,k] <- mean(a) 
    Upper[i,k] <- quantile(a, 0.95)
    Lower[i,k] <- quantile(a, 0.05)
    
    Mean_Traits[i,k] <- mean(b)
    Mean_Simp[i,k] <- mean(d)
    
  }
}

Mean <- apply(Mean[, -1], 1, mean)
Upper <- apply(Upper[, -1], 1, mean)
Lower <- apply(Lower[, -1], 1, mean)

Mean_Traits <- apply(Mean_Traits[, -1], 1, mean)
Mean_Simp <- apply(Mean_Simp[, -1], 1, mean)


# Random networks
MeanRand <-      Mean[which(seq$type == "random")]
UpperRand <-    Upper[which(seq$type == "random")]
LowerRand <-    Lower[which(seq$type == "random")]
Mean_TraitsRand <- Mean_Traits[which(seq$type == "random")]
Mean_DivRand <- Mean_Simp[which(seq$type == "random")]


#Scale free
MeanScale <- Mean[which(seq$type == "scalefree")]
UpperScale <- Upper[which(seq$type == "scalefree")]
LowerScale <- Lower[which(seq$type == "scalefree")]
Mean_TraitsScale <- Mean_Traits[which(seq$type == "scalefree")]
Mean_DivScale <- Mean_Simp[which(seq$type == "scalefree")]


# Small world

MeanSmall <- Mean[which(seq$type == "smallworld")]
UpperSmall <- Upper[which(seq$type == "smallworld")]
LowerSmall <- Lower[which(seq$type == "smallworld")]
Mean_TraitsSmall <- Mean_Traits[which(seq$type == "smallworld")]
Mean_DivSmall <- Mean_Simp[which(seq$type == "smallworld")]


  
  library(scales)
  #color stuff
  require(RColorBrewer)#load package
  col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values
  seqoverall <- seq
  
  
  #graphics.off()
  #png("NetworksWoDegrees.png", res = 900, height = 22, width = 20, units = "cm")
  
  par(mfrow=c(4,3),
      oma=c(2,2.6,1,0.5),
      mar=c(3,2,1.75,0.75))
  
  #Random network
  p <- 0.1
  g <- erdos.renyi.game(100, p, type = "gnp")
  V(g)$color <- "black"
  V(g)$frame.color <-"black"
  plot(g, vertex.label= NA, edge.arrow.size=0.1,vertex.size = 5)
  mtext("Random networks", side = 3, line = 0.5, cex = 1.1)
  legend("topright", "A", cex=1.1, bty="n")
  
  
  # Scale free networks
  
  g <- sample_pa(100, power = 1, m = 2, out.dist = NULL, out.seq = NULL,
                 out.pref = FALSE, zero.appeal = 1, directed = FALSE,
                 algorithm ="psumtree", start.graph = NULL)
  V(g)$color <- "black"
  V(g)$frame.color <- "black"
  plot(g, vertex.label= NA, edge.arrow.size=0.5,vertex.size = 5)
  mtext("Scale-free networks", side = 3, line = 0.5, cex = 1.1)
  legend("topright", "B", cex=1.1, bty="n")
  
  
  # Small world networks
  
  g <- watts.strogatz.game(1, 100, 4, p = 0.02, loops = FALSE, multiple = FALSE)
  V(g)$color <- "black"
  V(g)$frame.color <- "black"
  plot(g, vertex.label= NA, edge.arrow.size=0.5,vertex.size = 5)
  mtext("Small-world networks", side = 3, line = 0.5, cex = 1.1)
  legend("topright", "C", cex=1.1, bty="n")
  
  
  mu <- 1e-4
  
  seq <- seqoverall[which(seqoverall$type=="random"),]
  #Random network
  plot(MeanRand[which(seq$mu==mu)], ylim = c(0, 2000), xaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n", type = "n")
  par(xpd=FALSE)
  abline(h = 1000, lty = 2)
  axis(1,at=seq(1,10,1),labels=seq(0.1,1,0.1))
  arrows(seq(1,10,1),LowerRand[which(seq$mu==mu)],seq(1,10,1),UpperRand[which(seq$mu==mu)], code=3, lwd=2,col="grey", cex=1.3, length=0, angle = 90)
  points(MeanRand[which(seq$mu==mu)], pch=18, cex=2, col="black")
  mtext(expression(paste("Effective population size  ", italic(paste(N[e])))), side = 2, line = 2.8, cex = 0.9)
  legend("topright", "D", cex=1.1, bty="n")
  
  seq <- seqoverall[which(seqoverall$type=="scalefree"),]
  
  #Scale free network
  plot(MeanScale[which(seq$mu==mu)], ylim = c(0, 2000), type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n")
  par(xpd=FALSE)
  abline(h = 1000, lty = 2)
  axis(1,at=seq(1,10,1),labels=seq(0,1.5,length.out = 10))
  arrows(seq(1,10,1),LowerScale[which(seq$mu==mu)],seq(1,10,1),UpperScale[which(seq$mu==mu)], code=3, lwd=2,col="grey", cex=1.3, length=0, angle = 90)
  points(MeanScale[which(seq$mu==mu)], pch=18, cex=2, col="black")
  
  legend("topright", "E", cex=1.1, bty="n")
  
  seq <- seqoverall[which(seqoverall$type=="smallworld"),]
  
  #Small world network
  plot(MeanSmall[which(seq$mu==mu)], ylim = c(0, 2000), type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n")
  par(xpd=FALSE)
  abline(h = 1000, lty = 2)
  axis(1,at=seq(1,10,1),labels=seq(1,10,1))
  arrows(seq(1,10,1),LowerSmall[which(seq$mu==mu)],seq(1,10,1),UpperSmall[which(seq$mu==mu)], code=3, lwd=2,col="grey", cex=1.3, length=0, angle = 90)
  points(MeanSmall[which(seq$mu==mu)], pch=18, cex=2, col="black")
  
  legend("topright", "F", cex=1.1, bty="n")
  
  seq <- seqoverall[which(seqoverall$type=="random"),]
  
  
  #Random network
  plot(Mean_TraitsRand[which(seq$mu == mu)], log="y",ylim = c(1, 10000),type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
  axis(1,at=seq(1,10,1),labels=seq(0.1,1,0.1))
  axis(2,at=c(1,10,100,1000),labels=c(1,10,100,1000))
  
  for (mu in  unique(seq$mu)) {
    lines(Mean_TraitsRand[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
  }
  
  legend("topright", "G", cex=1.1, bty="n")
  
  mtext("Number of unique variants", side = 2, line = 3, cex = 0.9)
  legend("topleft",title = expression(paste("Innovation rate ", mu)), ncol = 4, legend=expression(10^-1,10^-2,10^-3,10^-4), col=c(col.pal[1:4]), pch = c(15,16,17,18), bty = "n", cex = 1.1)
  
  
  seq <- seqoverall[which(seqoverall$type=="scalefree"),]
  
  #Scale Free network
  plot(Mean_TraitsScale[which(seq$mu == mu)], log="y",ylim = c(1, 10000),type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
  axis(1,at=seq(1,10,1),labels=seq(0,1.5,length.out = 10))
  
  for (mu in  unique(seq$mu)) {
    lines(Mean_TraitsScale[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
  }
  legend("topright", "H", cex=1.1, bty="n")
  
  seq <- seqoverall[which(seqoverall$type=="smallworld"),]
  
  #Small world network
  plot(Mean_TraitsSmall[which(seq$mu == mu)], log="y", yaxt = "n",ylim = c(1, 10000),type = "n", xaxt = "n",  ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
  axis(1,at=seq(1,10,1),labels=seq(1,10,1))
  
  for (mu in  unique(seq$mu)) {
    lines(Mean_TraitsSmall[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
  }
  legend("topright", "I", cex=1.1, bty="n")
  
  
  
  #Random network
  seq <- seqoverall[which(seqoverall$type=="random"),]
  
  plot(Mean_DivRand[which(seq$mu == mu)],ylim = c(0, 1.2),type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
  axis(1,at=seq(1,10,1),labels=seq(0.1,1,0.1))
  axis(2,at=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1))
  
  for (mu in  unique(seq$mu)) {
    lines(Mean_DivRand[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
  }
  
  legend("topright", "J", cex=1.1, bty="n")
  mtext(expression(italic(p)), side = 1, line = 2.8, cex = 1.1)
  
  mtext("Simpson Diversity Index", side = 2, line = 3, cex = 0.9)
  
  
  seq <- seqoverall[which(seqoverall$type=="scalefree"),]
  
  #Scale Free network
  plot(Mean_DivScale[which(seq$mu == mu)],ylim = c(0, 1.2),type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
  axis(1,at=seq(1,10,1),labels=seq(0,1.5,length.out = 10))
  
  for (mu in  unique(seq$mu)) {
    lines(Mean_DivScale[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
  }
  legend("topright", "K", cex=1.1, bty="n")
  mtext(expression(italic(pi)), side = 1, line = 2.8, cex = 1.1)
  
  
  
  seq <- seqoverall[which(seqoverall$type=="smallworld"),]
  
  #Small world network
  plot(Mean_DivSmall[which(seq$mu == mu)], yaxt = "n",ylim = c(0, 1.2),type = "n", xaxt = "n",  ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
  axis(1,at=seq(1,10,1),labels=seq(1,10,1))
  
  for (mu in  unique(seq$mu)) {
    lines(Mean_DivSmall[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
  }
  legend("topright", "L", cex=1.1, bty="n")
  mtext(expression(italic(K)), side = 1, line = 2.8, cex = 1.1)
  
  
 # dev.off()
  

