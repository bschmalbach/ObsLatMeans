#observed score vs latent score, gender diffs
library(lavaan)
library(semTools)

df <- read.csv("data.csv",sep = "\t")
df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[,8:57]==0)==0,]

#recode items
rec_vec <- c(2,4,6,8,10,12,14,22,24,26,28:30,32,34,36,38,42,44,46)
for (i in rec_vec) {
  df[,i+7] <- (df[,i+7]-6)*(-1)
}

for (i in 1:50) {
  df[[i]] <- ordered(df[[i]], levels=sort(unique(df[[i]])))
}

m1 <- "O =~ O1 + O2 + O3 + O4 + O5 + O6 + O7 + O8 + O9 + O10"
m2 <- "C =~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10"
m3 <- "E =~ E1 + E2 + E3 + E4 + E5 + E6 + E7 + E8 + E9 + E10"
m4 <- "A =~ A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10"
m5 <- "N =~ N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9 + N10"

getLat <- function(x){
  fit <- sem(x, df, estimator="wlsmv", ordered=names(df)[8:57], group = "gender", group.equal=c("thresholds", "loadings", "intercepts"), std.lv=F, auto.fix.first=T)
  omega <- 0
  omega <- mean(compRelSEM(fit, ord.scale = T))
  if (is.na(omega)) {
    l1<-lavInspect(fit, what="std")$`1`$lambda; l2<-lavInspect(fit, what="std")$`2`$lambda
    omega <- mean(sum(abs(l1))^2 / (sum(abs(l1))^2 + sum(1-abs(l1)^2)), sum(abs(l2))^2 / (sum(abs(l2))^2 + sum(1-abs(l2^2))))
  }
  fitM <- fitmeasures(fit, fit.measures = c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  return(list(lavInspect(fit, what="std")$`2`$alpha, 
              omega,
              fitM))
}

lv <- list(getLat(m1), getLat(m2), getLat(m3), getLat(m4), getLat(m5))
lat <- -1*unlist(lapply(lv, `[[`, 1))

#obs diff
for (i in 1:50) {
  df[[i]] <- as.numeric(paste(df[[i]]))
}
df$O_sum <- scale(rowSums(df[c("O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8", "O9", "O10")]))
df$C_sum <- scale(rowSums(df[c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10")]))
df$E_sum <- scale(rowSums(df[c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10")]))
df$A_sum <- scale(rowSums(df[c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10")]))
df$N_sum <- scale(rowSums(df[c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10")]))

get_d <- function(x){
  v1 <- x[df$gender==1] 
  v2 <- x[df$gender==2]
  mean_diff <- mean(v1)-mean(v2)
  pool_sd <- sqrt((var(v1)+var(v2))/2)
  d <- mean_diff/pool_sd
  return(d)
}

obs <- c(get_d(df$O_sum),
         get_d(df$C_sum),
         get_d(df$E_sum),
         get_d(df$A_sum),
         get_d(df$N_sum))


res <- data.frame(cbind(lat, obs))
names(res)[1] <- "lat"

#compare

n_g <- c(sum(df$gender==1), sum(df$gender==2))

res$diff <- res$lat-res$obs
res$ratio <- res$obs/res$lat

res$var1 <- ((n_g[1]+n_g[2])/(n_g[1]*n_g[2])) + (res$lat^2 / (2*(n_g[1]+n_g[2])))
res$var2 <- ((n_g[1]+n_g[2])/(n_g[1]*n_g[2])) + (res$obs^2 / (2*(n_g[1]+n_g[2])))
res$pool_sd <- sqrt((res$var1+res$var2)/2)

res$z <- res$diff/(res$pool_sd)
res$p <- 2*(1-pnorm(abs(res$z)))
res$dd <- res$diff/(res$pool_sd)

res$UL_CI <- res$diff + 1.96*res$pool_sd
res$LL_CI <- res$diff - 1.96*res$pool_sd
res$n <- sum(n_g)
res$omega <- lapply(lv, `[[`, 2)
res$fitM <- lapply(lv, `[[`, 3)
