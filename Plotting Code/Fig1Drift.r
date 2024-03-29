
####
###
##
# Plotting code for Fig. 1
##
###
####


N_pop <- 8
 t_max <- 100
 p0 <- 0.5
 
 
 library(scales)
 #color stuff
 require(RColorBrewer)#load package
 col.pal <- brewer.pal(8, "Dark2") #create a pallette which you loop over for corresponding values
 seqoverall <- seq
 
 
 require(RColorBrewer)#load package
 col.pal <- brewer.pal(N_pop, "Dark1") #create a pallette which you loop over for corresponding values
 
 # graphics.off()
 # png("Drift.png", width = 15,height = 15, units = "cm", res = 1200)
 
 par(mfrow= c(2,2), 
     mar=c(1.4,2,2.2,0.4), 
     oma=c(2,2.2,0,0))
 
 N<- 10
 
 p <- matrix(0, N_pop, t_max)
 
 for (i in 1:N_pop) {
   Pop <- c(rep(1, p0*N), rep(0, (1-p0)*N))
   
   for (t in 1:t_max) {
     
     p[i,t] <- sum(Pop)/length(Pop)
     Pop <- sample(Pop, replace = TRUE)
   }
 }
 
 plot(p[1,], type = "n", ylim = c(0,1), ylab = "")
 for (i in 1:N_pop) {
 lines(p[i,], col = col.pal[i], lwd=2)
 }
 mtext(side = 3, expression(paste(italic("N")," = ",10)))
 
  N<- 100

  p <- matrix(0, N_pop, t_max)
  
  for (i in 1:N_pop) {
    Pop <- c(rep(1, p0*N), rep(0, (1-p0)*N))
    
    for (t in 1:t_max) {
      
      p[i,t] <- sum(Pop)/length(Pop)
      Pop <- sample(Pop, replace = TRUE)
    }
  }
  
  plot(p[1,], type = "n", ylim = c(0,1), ylab = "")
  for (i in 1:N_pop) {
    lines(p[i,], col = col.pal[i], lwd=2)
  }
  mtext(side = 3, expression(paste(italic("N")," = ",100)))
  
   N<- 1000
   p <- matrix(0, N_pop, t_max)
   
   for (i in 1:N_pop) {
     Pop <- c(rep(1, p0*N), rep(0, (1-p0)*N))
     
     for (t in 1:t_max) {
       
       p[i,t] <- sum(Pop)/length(Pop)
       Pop <- sample(Pop, replace = TRUE)
     }
   }
   
   plot(p[1,], type = "n", ylim = c(0,1), ylab = "")
   for (i in 1:N_pop) {
     lines(p[i,], col = col.pal[i], lwd=2)
   }
   mtext(side = 3, expression(paste(italic("N")," = ",1000)))
   
   
   N<- 10000
   
   p <- matrix(0, N_pop, t_max)
   
   for (i in 1:N_pop) {
      Pop <- c(rep(1, p0*N), rep(0, (1-p0)*N))
      
      for (t in 1:t_max) {
         
         p[i,t] <- sum(Pop)/length(Pop)
         Pop <- sample(Pop, replace = TRUE)
      }
   }
   
   plot(p[1,], type = "n", ylim = c(0,1), ylab = "")
   for (i in 1:N_pop) {
      lines(p[i,], col = col.pal[i], lwd=2)
   }
   mtext(side = 3, expression(paste(italic("N")," = ",10000)))
   mtext(side = 1, line = 1 , "Generation", outer = TRUE, cex = 1)
   mtext(side = 2, line = 1 , "Allele frequency", outer = TRUE, cex = 1)
   
   #dev.off()