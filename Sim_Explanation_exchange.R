
y <- list()
z1 <- list()
z2 <- list()

y[[1]] <- 1:10000
z1[[1]] <- 1:5000
z2[[1]] <- 5001:10000

Div_full <- c()
Div_sep <- c()

Div_full[1] <- length(unique(y[[1]]))
Div_sep[1] <- length(unique(c(z1[[1]])  )   )


x <- 10000

for (t in 2:50) {
  
y[[t]] <- sample(y[[t-1]], size = 10000, replace = TRUE)  

Div_full[t] <-  length(unique(y[[t]]))





z1[[t]] <- sample(c(z1[[t-1]],z2[[t-1]][1:5000]), size = 5000, replace = TRUE)  
z2[[t]] <- sample(c(z1[[t-1]][1:5000],z2[[t-1]]), size = 5000, replace = TRUE)  

Div_sep[t] <- length(unique(c(z1[[t]])  )   )



}

plot(Div_full, type = "l")
lines(Div_sep, col = "blue")
