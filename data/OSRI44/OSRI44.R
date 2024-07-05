#observed score vs latent score, gender diffs

library(lavaan)
set.seed(1)
df <- read.csv("data.csv",sep = "\t")
df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[1:44]==0)==0,]

#recode items
rec_vec <- c(9, 13, 15)
for (i in rec_vec) {
  df[,i] <- (6-df[,i])
}

model <- "MASC =~ Q1 + Q3 + Q5 + Q7 + Q9 + Q11 + Q13 + Q15 + Q17 + Q19 + Q21 + Q23 + Q25 + Q27 + Q29 + Q31 + Q33 + Q35 + Q37 + Q39 + Q41 + Q43
          FEM =~ Q2 + Q4 + Q6 + Q8 + Q10 + Q12 + Q14 + Q16 + Q18 + Q20 + Q22 + Q24 + Q26 + Q28 + Q30 + Q32 + Q34 + Q36 + Q38 + Q40 + Q42 + Q44"

fit <- sem(model, df, estimator="mlr", group = "gender", group.equal=c("loadings", "intercepts"), std.lv=T, auto.fix.first=F)
lat <- (lavInspect(fit, what="std")$`1`$alpha)


#obs diff
df$MASK_sum <- scale(rowSums(df[seq(1, 43, 2)]))
df$FEM_sum <- scale(rowSums(df[seq(2, 44, 2)]))

get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- c(get_d(df$MASK_sum),
         get_d(df$FEM_sum))

res <- data.frame(cbind(lat, obs))
names(res)[1] <- "lat"

#compare

n_g <- as.numeric(c(sum(df$gender==1), sum(df$gender==2)))

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
res$omega <- colMeans(matrix(unlist(semTools::reliability(fit, what="omega")),nrow = nrow(res), byrow = T))