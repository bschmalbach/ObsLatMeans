#observed score vs latent score, gender diffs
library(lavaan)

df <- read.csv("data.csv",sep = ",")
df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[1:32]==-1)==0,]

#recode items
rec_vec <- c(1,9,17,25,29,22,7,15,23,31,16)
for (i in rec_vec) {
  df[,i] <- (6-df[,i])
}

model <- "HSQ_affiliative =~ Q1 + Q5 + Q9 + Q13 + Q17 + Q21 + Q25 + Q29
          HSQ_selfenhancing =~ Q2 + Q6 + Q10 + Q14 + Q18 + Q22 + Q26 + Q30
          HSQ_aggressive =~ Q3 + Q7 + Q11 + Q15 + Q19 + Q23 + Q27 + Q31
          HSQ_selfdefeating =~ Q4 + Q8 + Q12 + Q16 + Q20 + Q24 + Q28 + Q32"

fit <- sem(model, df, estimator="mlr", group = "gender", group.equal=c("loadings", "intercepts"), std.lv=T, auto.fix.first=F)
lat <- (lavInspect(fit, what="std")$`1`$alpha)


#obs diff
df$affiliative_sum <- scale(rowSums(df[c("Q1", "Q5", "Q9", "Q13", "Q17", "Q21", "Q25", "Q29")]))
df$selfenhancing_sum <- scale(rowSums(df[c("Q2", "Q6", "Q10", "Q14", "Q18", "Q22", "Q26", "Q30")]))
df$aggressive_sum <- scale(rowSums(df[c("Q3", "Q7", "Q11", "Q15", "Q19", "Q23", "Q27", "Q31")]))
df$selfdefeating_sum <- scale(rowSums(df[c("Q4", "Q8", "Q12", "Q16", "Q20", "Q24", "Q28", "Q32")]))

get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- c(get_d(df$affiliative_sum),
         get_d(df$selfenhancing_sum),
         get_d(df$aggressive_sum),
         get_d(df$selfdefeating_sum))

res <- data.frame(cbind(lat, obs))
names(res)[1] <- "lat"

#compare

n_g <- c(sum(df$gender==1), sum(df$gender==2))

res$diff <- res$lat-res$obs
res$ratio <- res$obs/res$lat

res$var1 <- ((n_g[1]+n_g[2])/(n_g[1]*n_g[2])) + (res$lat^2 / (2*(n_g[1]+n_g[2])))
res$var2 <- ((n_g[1]+n_g[2])/(n_g[1]*n_g[2])) + (res$obs^2 / (2*(n_g[1]+n_g[2])))
res$pool_sd <- sqrt((res$var1+res$var2)/2)

res$z <- res$diff/(res$pool_sd/sqrt(sum(n_g)))
res$p <- 2*(1-pnorm(abs(res$z)))
res$dd <- res$diff/(res$pool_sd)

res$UL_CI <- res$diff + 1.96*res$pool_sd
res$LL_CI <- res$diff - 1.96*res$pool_sd

res$n <- sum(n_g)
res$omega <-  rowMeans(matrix(unlist(semTools::reliability(fit, what="omega")),nrow = nrow(res), byrow = F))