
# Effective population size analysis script
# Social network structure

library(igraph)

# Random networks

setwd("C:/Users/dominik_deffner/Documents/")
load("~/CulturalEffectivePopulationSize/Neff_randnet_0909")
seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-4), p = seq(0.1,1,0.1) )

MeanRand <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
UpperRand <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
LowerRand <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
    }
    MeanRand[i,k] <- mean(a) 
    UpperRand[i,k] <- quantile(a, 0.95)
    LowerRand[i,k] <- quantile(a, 0.05)
    
  }
}

MeanRand <- apply(MeanRand[, -1], 1, mean)
UpperRand <- apply(UpperRand[, -1], 1, mean)
LowerRand <- apply(LowerRand[, -1], 1, mean)



#Scale free
load("~/CulturalEffectivePopulationSize/Neff_scalefreenet_0909")
seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-4), p = seq(0,2,0.1) )


MeanScale <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
UpperScale <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])
LowerScale <- matrix(NA, nrow = nrow(seq), ncol = seq$tmax[1])

for (i in 1 : nrow(seq)) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
    }
    MeanScale[i,k] <- mean(a) 
    UpperScale[i,k] <- quantile(a, 0.95)
    LowerScale[i,k] <- quantile(a, 0.05)
    
  }
}

MeanScale <- apply(MeanScale[, -1], 1, mean)
UpperScale <- apply(UpperScale[, -1], 1, mean)
LowerScale <- apply(LowerScale[, -1], 1, mean)



# Small world

load("~/CulturalEffectivePopulationSize/Neff_smallworldnet_0909")
seq<-expand.grid(N=1000, tmax=300,Nsim = 1000, mu = c(1e-4), p = c(0,0.01,0.1,0.5), m = seq(2,20,2) )

xx <- which(seq$p == 0.01)

MeanSmall <- matrix(NA, nrow = length(xx), ncol = seq$tmax[1])
UpperSmall <- matrix(NA, nrow = length(xx), ncol = seq$tmax[1])
LowerSmall <- matrix(NA, nrow = length(xx), ncol = seq$tmax[1])

for (i in xx) {
  print(i)
  for (k in 1: seq$tmax[1]) {
    a <- c()
    for (j in 1:  seq$Nsim[1]) {
      a <- c(a, result[[i]][[j]][["N_effective"]][[1]][k])
    }
    MeanSmall[which(xx == i),k] <- mean(a) 
    UpperSmall[which(xx == i),k] <- quantile(a, 0.95)
    LowerSmall[which(xx == i),k] <- quantile(a, 0.05)
    
  }
}

MeanSmall <- apply(MeanSmall[, -1], 1, mean)
UpperSmall <- apply(UpperSmall[, -1], 1, mean)
LowerSmall <- apply(LowerSmall[, -1], 1, mean)












library(scales)
#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values


graphics.off()
png("Networks.png", res = 1200, height = 20, width = 20, units = "cm")

par(mfrow=c(3,3),
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



#Random network
plot(MeanRand, ylim = c(0, 1500), xaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n", type = "n")
par(xpd=FALSE)
abline(h = 1000, lty = 2)
axis(1,at=seq(1,10,1),labels=seq(0.1,1,0.1))
arrows(seq(1,10,1),LowerRand,seq(1,10,1),UpperRand, code=3, lwd=2,col=col.pal[1], cex=1.3, length=0, angle = 90)
points(MeanRand, pch=18, cex=2, col=col.pal[1])
mtext(expression(italic(p)), side = 1, line = 2.8, cex = 1.1)
mtext(expression(paste("Effective population size  ", italic(paste(N[e])))), side = 2, line = 3, cex = 1.1)
legend("topright", "G", cex=1.1, bty="n")

#Scale free network
plot(MeanScale, ylim = c(0, 1500), type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n")
par(xpd=FALSE)
abline(h = 1000, lty = 2)
axis(1,at=seq(1,21,2),labels=seq(0,2,0.2))
arrows(seq(1,21,1),LowerScale,seq(1,21,1),UpperScale, code=3, lwd=2,col=col.pal[2], cex=1.3, length=0, angle = 90)
points(MeanScale, pch=18, cex=2, col=col.pal[2])
mtext(expression(italic(pi)), side = 1, line = 2.8, cex = 1.1)

legend("topright", "H", cex=1.1, bty="n")


#Small world network
plot(MeanSmall, ylim = c(0, 1500), type = "n", xaxt = "n", yaxt = "n", ylab = "", pch=18, cex=1.3, xlab = "", main = "", bty="n")
par(xpd=FALSE)
abline(h = 1000, lty = 2)
axis(1,at=seq(1,10,1),labels=seq(2,20,2))
arrows(seq(1,21,1),LowerSmall,seq(1,21,1),UpperSmall, code=3, lwd=2,col=col.pal[3], cex=1.3, length=0, angle = 90)
points(MeanSmall, pch=18, cex=2, col=col.pal[3])
mtext(expression(italic(K)), side = 1, line = 2.8, cex = 1.1)

legend("topright", "I", cex=1.1, bty="n")

dev.off()
