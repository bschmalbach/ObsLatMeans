#observed score vs latent score, gender diffs

library(lavaan)

df <- read.csv("data.csv",sep = "\t")
df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[11:22]==0)==0,]

for (i in 11:22) {
  df[[i]] <- ordered(df[[i]], levels=sort(unique(df[[i]])))
}

m1 <- "DDP =~ DDP1 + DDP2 + DDP3 + DDP4"
m2 <- "DDN =~ DDN1 + DDN2 + DDN3 + DDN4"
m3 <- "DDM =~ DDM1 + DDM2 + DDM3 + DDM4"

#lat diff
getLat <- function(x){
  fit <- sem(x, df, estimator="wlsmv", ordered=names(df)[11:22], group = "gender", group.equal=c("thresholds", "loadings", "intercepts"), std.lv=F, auto.fix.first=T)
  omega <- 0
  omega <- mean(semTools::compRelSEM(fit, ord.scale = T))
  if (is.na(omega)) {
    l1<-lavInspect(fit, what="std")$`1`$lambda; l2<-lavInspect(fit, what="std")$`2`$lambda
    omega <- mean(sum(abs(l1))^2 / (sum(abs(l1))^2 + sum(1-abs(l1)^2)), sum(abs(l2))^2 / (sum(abs(l2))^2 + sum(1-abs(l2^2))))
  }
  fitM <- fitmeasures(fit, fit.measures = c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  return(list(lavInspect(fit, what="std")$`1`$alpha, 
              omega,
              fitM))
}

lv <- list(getLat(m1), getLat(m2), getLat(m3))
lat <- unlist(lapply(lv, `[[`, 1))


#obs diff
for (i in 11:22) {
  df[[i]] <- as.numeric(paste(df[[i]]))
}
df$DDP_sum <- scale(rowSums(df[11:14]))
df$DDN_sum <- scale(rowSums(df[15:18]))
df$DDM_sum <- scale(rowSums(df[19:22]))

get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- c(get_d(df$DDP_sum),
         get_d(df$DDN_sum),
         get_d(df$DDM_sum))

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