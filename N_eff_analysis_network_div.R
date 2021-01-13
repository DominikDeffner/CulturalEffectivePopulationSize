
# Effective population size analysis script
# Social network structure

library(igraph)

# Random networks

load("/home/dominik_deffner/Neff_randnet_1612")
seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-1,1e-2,1e-3,1e-4), p = seq(0.1,1,0.1) )

MeanRand <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
UpperRand <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
LowerRand <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_TraitsRand <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    b <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
      b <- c(b, result[[i]][[j]][["Div_NoTraits"]][[1]][k])
      
    }
    MeanRand[i,k] <- mean(a) 
    UpperRand[i,k] <- quantile(a, 0.95)
    LowerRand[i,k] <- quantile(a, 0.05)
    
    Mean_TraitsRand[i,k] <- mean(b)
    
  }
}

MeanRand <- apply(MeanRand[, -1], 1, mean)
UpperRand <- apply(UpperRand[, -1], 1, mean)
LowerRand <- apply(LowerRand[, -1], 1, mean)

Mean_TraitsRand <- apply(Mean_TraitsRand[, -1], 1, mean)



#Scale free
load("/home/dominik_deffner/Neff_scalefree_1612")
seq<-expand.grid(N=1000, tmax=100,Nsim = 50, mu = c(1e-1,1e-2,1e-3,1e-4), pi = seq(0,1.5,length.out = 10) )


MeanScale <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
UpperScale <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
LowerScale <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_TraitsScale <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    b <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
      b <- c(b, result[[i]][[j]][["Div_NoTraits"]][[1]][k])
      
    }
    MeanScale[i,k] <- mean(a) 
    UpperScale[i,k] <- quantile(a, 0.95)
    LowerScale[i,k] <- quantile(a, 0.05)
    
    Mean_TraitsScale[i,k] <- mean(b)
  }
}

MeanScale <- apply(MeanScale[, -1], 1, mean)
UpperScale <- apply(UpperScale[, -1], 1, mean)
LowerScale <- apply(LowerScale[, -1], 1, mean)

Mean_TraitsScale <- apply(Mean_TraitsScale[, -1], 1, mean)



# Small world

load("/home/dominik_deffner/Neff_smallworldnet_1612")
seq<-expand.grid(N=1000, tmax=100,Nsim = 50, mu = c(1e-1,1e-2,1e-3,1e-4), K = seq(1,10,1) )


MeanSmall <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
UpperSmall <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
LowerSmall <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

Mean_TraitsSmall <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    b <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
      b <- c(b, result[[i]][[j]][["Div_NoTraits"]][[1]][k])
      
    }
    MeanSmall[i,k] <- mean(a) 
    UpperSmall[i,k] <- quantile(a, 0.95)
    LowerSmall[i,k] <- quantile(a, 0.05)
    
    Mean_TraitsSmall[i,k] <- mean(b)
  }
}

MeanSmall <- apply(MeanSmall[, -1], 1, mean)
UpperSmall <- apply(UpperSmall[, -1], 1, mean)
LowerSmall <- apply(LowerSmall[, -1], 1, mean)

Mean_TraitsSmall <- apply(Mean_TraitsSmall[, -1], 1, mean)












library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values


graphics.off()
png("Networks.png", res = 400, height = 26, width = 20, units = "cm")

par(mfrow=c(4,3),
    oma=c(2,2.6,1,0.5),
    mar=c(3,2,1.75,0.75))

#Random network
p <- 0.1
g <- erdos.renyi.game(100, p, type = "gnp")
V(g)$color <- col.pal[1]
V(g)$frame.color <- col.pal[1]
plot(g, vertex.label= NA, edge.arrow.size=0.1,vertex.size = 5)
mtext("Random networks", side = 3, line = 0.5, cex = 1.1)
legend("topright", "A", cex=1.1, bty="n")



# Scale free networks

g <- sample_pa(100, power = 1, m = 1, out.dist = NULL, out.seq = NULL,
               out.pref = FALSE, zero.appeal = 1, directed = FALSE,
               algorithm ="psumtree", start.graph = NULL)
V(g)$color <- col.pal[2]
V(g)$frame.color <- col.pal[2]
plot(g, vertex.label= NA, edge.arrow.size=0.5,vertex.size = 5)
mtext("Scale-free networks", side = 3, line = 0.5, cex = 1.1)
legend("topright", "B", cex=1.1, bty="n")


# Small world networks

g <- watts.strogatz.game(1, 100, 4, p = 0.02, loops = FALSE, multiple = FALSE)
V(g)$color <- col.pal[3]
V(g)$frame.color <- col.pal[3]
plot(g, vertex.label= NA, edge.arrow.size=0.5,vertex.size = 5)
mtext("Small-world networks", side = 3, line = 0.5, cex = 1.1)
legend("topright", "C", cex=1.1, bty="n")



# Degree distribution
#Random
aa <- list()
for (i in 1:1000) {
  g <- erdos.renyi.game(1000, 0.1, type = "gnp")
  aa[[i]] <- degree.distribution(g)
}

plot("", xlab = "", ylab = "", log = "xy", type = "n", xlim = c(1e-0,1e3), ylim = c(1e-3,1e-0), bty = "n", xaxt = "n", yaxt ="n")
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
axis(2,at=c(1e-3,1e-2, 1e-1, 1e-0),labels=c(1e-3,1e-2, 1e-1, 1e-0))
legend("topright", "D", cex=1.1, bty="n")

for (i in 1:10000) {
  lines(aa[[i]], log = "xy", col= col.pal[1], pch = 16, lwd=1.3)
}

mtext("P(k)", side = 2, line = 3, cex = 1.1)


#Scale - Free
aa <- list()
for (i in 1:1000) {
  g <- sample_pa(1000, power = 1, m = 1, out.dist = NULL, out.seq = NULL,out.pref = FALSE, zero.appeal = 1, directed = FALSE,algorithm ="psumtree", start.graph = NULL) 
  aa[[i]] <- degree.distribution(g)
}

plot("", xlab = "", ylab = "", log = "xy", type = "n", xlim = c(1e-0,1e3), ylim = c(1e-3,1e-0), bty = "n", xaxt = "n", yaxt ="n")
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))
legend("topright", "E", cex=1.1, bty="n")

for (i in 1:1000) {
  lines(aa[[i]], log = "xy", col= col.pal[2], pch = 16, lwd=1.3)
}

mtext(paste("Degree ","k"), side = 1, line = 3, cex = 1.1)

# Small world
aa <- list()
for (i in 1:1000) {
  g <- watts.strogatz.game(1, 1000, 4, p = 0.02, loops = FALSE, multiple = FALSE)
  aa[[i]] <- degree.distribution(g)
}

plot("", xlab = "", ylab = "", log = "xy", type = "n", xlim = c(1e-0,1e3), ylim = c(1e-3,1e-0), bty = "n", xaxt = "n", yaxt ="n")
axis(1,at=c(1,10,100,1000),labels=c(1,10,100,1000))

for (i in 1:1000) {
  lines(aa[[i]], log = "xy", col= col.pal[3], pch = 16, lwd=1.3)
}
legend("topright", "F", cex=1.1, bty="n")


mu <- 1e-4

#Random network
plot(MeanRand[which(seq$mu==mu)], ylim = c(0, 2000), xaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n", type = "n")
par(xpd=FALSE)
abline(h = 1000, lty = 2)
axis(1,at=seq(1,10,1),labels=seq(0.1,1,0.1))
arrows(seq(1,10,1),LowerRand[which(seq$mu==mu)],seq(1,10,1),UpperRand[which(seq$mu==mu)], code=3, lwd=2,col=col.pal[1], cex=1.3, length=0, angle = 90)
points(MeanRand, pch=18, cex=2, col=col.pal[1])
mtext(expression(italic(p)), side = 1, line = 2.8, cex = 1.1)
mtext(expression(paste("Effective population size  ", italic(paste(N[e])))), side = 2, line = 3, cex = 1.1)
legend("topright", "G", cex=1.1, bty="n")

#Scale free network
plot(MeanScale[which(seq$mu==mu)], ylim = c(0, 2000), type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n")
par(xpd=FALSE)
abline(h = 1000, lty = 2)
axis(1,at=seq(1,10,1),labels=seq(0,1.5,length.out = 10))
arrows(seq(1,10,1),LowerScale[which(seq$mu==mu)],seq(1,10,1),UpperScale[which(seq$mu==mu)], code=3, lwd=2,col=col.pal[2], cex=1.3, length=0, angle = 90)
points(MeanScale[which(seq$mu==mu)], pch=18, cex=2, col=col.pal[2])
mtext(expression(italic(pi)), side = 1, line = 2.8, cex = 1.1)

legend("topright", "H", cex=1.1, bty="n")


#Small world network
plot(MeanSmall[which(seq$mu==mu)], ylim = c(0, 2000), type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n")
par(xpd=FALSE)
abline(h = 1000, lty = 2)
axis(1,at=seq(1,10,1),labels=seq(1,10,1))
arrows(seq(1,10,1),LowerSmall[which(seq$mu==mu)],seq(1,10,1),UpperSmall[which(seq$mu==mu)], code=3, lwd=2,col=col.pal[3], cex=1.3, length=0, angle = 90)
points(MeanSmall[which(seq$mu==mu)], pch=18, cex=2, col=col.pal[3])
mtext(expression(italic(K)), side = 1, line = 2.8, cex = 1.1)

legend("topright", "I", cex=1.1, bty="n")


#Random network
plot(Mean_TraitsRand[which(seq$mu == mu)], log="y",ylim = c(1, 10000),type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
axis(1,at=seq(1,10,1),labels=seq(0.1,1,0.1))
axis(2,at=c(1,10,100,1000),labels=c(1,10,100,1000))

for (mu in  unique(seq$mu)) {
  lines(Mean_TraitsRand[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
}

legend("topright", "J", cex=1.1, bty="n")
mtext(expression(italic(p)), side = 1, line = 2.8, cex = 1.1)

mtext("Number of unique variants", side = 2, line = 3, cex = 1.1)



#Scale Free network
plot(Mean_TraitsScale[which(seq$mu == mu)], log="y",ylim = c(1, 10000),type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
axis(1,at=seq(1,10,1),labels=seq(0,1.5,length.out = 10))

for (mu in  unique(seq$mu)) {
  lines(Mean_TraitsScale[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
}
legend("topright", "K", cex=1.1, bty="n")
mtext(expression(italic(pi)), side = 1, line = 2.8, cex = 1.1)


#Small world network
plot(Mean_TraitsSmall[which(seq$mu == mu)], log="y", yaxt = "n",ylim = c(1, 10000),type = "n", xaxt = "n",  ylab = "", pch=18, cex=1.1, xlab = "", bty = "n")
axis(1,at=seq(1,10,1),labels=seq(1,10,1))

for (mu in  unique(seq$mu)) {
  lines(Mean_TraitsSmall[which(seq$mu == mu)],type = "b",  pch=14+which(unique(seq$mu)==mu), cex=1.2, col= col.pal[which(unique(seq$mu)==mu)])
}
legend("topright", "L", cex=1.1, bty="n")
mtext(expression(italic(K)), side = 1, line = 2.8, cex = 1.1)

legend("topright",title = expression(paste("Innovation rate ", mu)), ncol = 2, legend=expression(10^-1,10^-2,10^-3,10^-4), col=c(col.pal[1:4]), pch = c(15,16,17,18), bty = "n", cex = 1.3)



dev.off()
