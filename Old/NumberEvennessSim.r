
# Just trying to get some understanding of N_eff

N_e <- function(k_bar,V_k){
  (N * k_bar - 1 ) / ( (V_k / k_bar) + k_bar  -1 )
}

# k_bar = 1
N_eff <- function(V_k){
  (N * 1 - 1 ) / ( (V_k / 1) + 1  -1 )
}



library(rethinking)
library(vegan)

N <- 1000
N_variants <- 1000
Div <- seq(1/N_variants,1,0.01)


par (mfrow = c(1,5), 
     oma=c(3,2.5,0,0),
     mar=c(1.5,2.4,1.5,0.1))


for (f in c(0.1, 0.5, 1, 1.5, 3)) {
  


result_Ne <- c()
result_Div <- c()

for (i in Div) {
print(which(Div == i))
    
p <- c()
p[1] <- i
p[2:N_variants] <- (1-i)/(N_variants-1)

 
N_sim <- 100
Ne <- c()
Diversity <- c()

for (sim in 1:N_sim) {

Pop <- sample(1:N_variants, N, replace = TRUE, prob = p)

Diversity[sim] <- diversity(sapply(unique(Pop), function (x) length(which(Pop == x))), index = "simpson")


#Frequency of each variant
Variants <- unique(Pop)

Freq_Variants <- sapply(Variants, function (x) length(which(Pop == x)))

#Probability individuals choose each variant
P <- Freq_Variants^f / sum(Freq_Variants^f)
P_Ind <- P/Freq_Variants
Copied <- sample(1:N, N, replace = TRUE, sapply(1:N, function (x) P_Ind[which(Variants == Pop[x])]))

Offspring_Record<- sapply(1:N, function (x) length(which(Copied == x)))

Ne[sim] <- N_eff(var(Offspring_Record))

}

result_Ne[which(Div == i)] <- mean(Ne)
result_Div[which(Div == i)] <- mean(Diversity)

}

#plot(Div, result_Ne, type = "b", ylim = c(0, N + 0.1*N), main = paste("Theta" = f)) #expression(paste("theta = ", f)))
plot(result_Div, result_Ne, type = "b", ylim = c(0, N + 0.1*N), main = paste("Theta" = f)) #expression(paste("theta = ", f)))

abline(h = N, lty = 2)

}
mtext("Diversity", line = 1.3, side = 1, outer = TRUE)

#mtext("Frequency of 1st variant", line = 1.3, side = 1, outer = TRUE)
mtext("Effective Populatio Size", line = 1.3, side = 2, outer = TRUE)




