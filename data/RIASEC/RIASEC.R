#observed score vs latent score, gender diffs
library(lavaan)

df <- read.csv("data.csv",sep = "\t")
df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[1:48]==0)==0,]
for (i in 1:48) {
  df[[i]] <- ordered(df[[i]], levels=sort(unique(df[[i]])))
}

m1 <- "R =~ R1 + R2 + R3 + R4 + R5 + R6 + R7 + R8"
m2 <- "I =~ I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8"
m3 <- "A =~ A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8"
m4 <- "S =~ S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8"
m5 <- "E =~ E1 + E2 + E3 + E4 + E5 + E6 + E7 + E8"
m6 <- "C =~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8"

getLat <- function(x){
  fit <- sem(x, df, estimator="wlsmv", ordered=names(df)[1:48], group = "gender", group.equal=c("thresholds", "loadings", "intercepts"), std.lv=F, auto.fix.first=T)
  omega <- 0
  omega <- mean(semTools::compRelSEM(fit, ord.scale = T))
  if (is.na(omega)) {
    l1<-lavInspect(fit, what="std")$`1`$lambda; l2<-lavInspect(fit, what="std")$`2`$lambda
    omega <- mean(sum(abs(l1))^2 / (sum(abs(l1))^2 + sum(1-abs(l1)^2)), sum(abs(l2))^2 / (sum(abs(l2))^2 + sum(1-abs(l2^2))))
  }
  fitM <- fitmeasures(fit, fit.measures = c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  return(list(lavInspect(fit, what="std")$`2`$alpha*-1, 
              omega,
              fitM))
}

lv <- list(getLat(m1), getLat(m2), getLat(m3), getLat(m4), getLat(m5), getLat(m6))
lat <- unlist(lapply(lv, `[[`, 1))
#obs diff
for (i in 1:48) {
  df[[i]] <- as.numeric(paste(df[[i]]))
}
df$R_sum <- scale(rowSums(df[1:8]))
df$I_sum <- scale(rowSums(df[9:16]))
df$A_sum <- scale(rowSums(df[17:24]))
df$S_sum <- scale(rowSums(df[25:32]))
df$E_sum <- scale(rowSums(df[33:40]))
df$C_sum <- scale(rowSums(df[41:48]))



get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- c(get_d(df$R_sum),
         get_d(df$I_sum),
         get_d(df$A_sum),
         get_d(df$S_sum),
         get_d(df$E_sum),
         get_d(df$C_sum))

res <- data.frame(cbind(lat, obs))
names(res)[1] <- "lat"

#compare
n_g <- c(sum(df$gender==1), sum(df$gender==2))

res$diff <- res$lat-res$obs
res$ratio <- res$obs/res$lat

res$var1 <- ((n_g[1]+n_g[2])/(1.0*n_g[1]*n_g[2])) + (res$lat^2 / (2*(n_g[1]+n_g[2])))
res$var2 <- ((n_g[1]+n_g[2])/(1.0*n_g[1]*n_g[2])) + (res$obs^2 / (2*(n_g[1]+n_g[2])))
res$pool_sd <- sqrt((res$var1+res$var2)/2)

res$z <- res$diff/(res$pool_sd)
res$p <- 2*(1-pnorm(abs(res$z)))
res$dd <- res$diff/(res$pool_sd)

res$UL_CI <- res$diff + 1.96*res$pool_sd
res$LL_CI <- res$diff - 1.96*res$pool_sd
res$n <- sum(n_g)
res$omega <- lapply(lv, `[[`, 2)
res$fitM <- lapply(lv, `[[`, 3)