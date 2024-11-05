#observed score vs latent score, gender diffs
library(lavaan)

df <- read.csv("data.csv")
df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[,1:40]==0)==0,]

#recode items
model <- "AS =~ AS1 + AS2 + AS3 + AS4 + AS5 + AS6 + AS7 + AS8 + AS9 + AS10
  SC =~ SC1 + SC2 + SC3 + SC4 + SC5 + SC6 + SC7 + SC8 + SC9 + SC10
  AD =~ AD1 + AD2 + AD3 + AD4 + AD5 + AD6 + AD7 + AD8 + AD9 + AD10
  DO =~ DO1 + DO2 + DO3 + DO4 + DO5 + DO6 + DO7 + DO8 + DO9 + DO10"
fit <- sem(model, df, estimator="wlsmv", ordered=names(df)[1:40])
rowSums(lavInspect(fit, what="std")$lambda)
rec_vec <- 1:40
rec_vec <- rec_vec[-rowSums(lavInspect(fit, what="std")$lambda)>0]
for (i in rec_vec) {
  df[,i] <- 6-df[,i]
}

for (i in 1:40) {
  df[[i]] <- ordered(df[[i]], levels=sort(unique(df[[i]])))
}

m1 <- "AS =~ AS1 + AS2 + AS3 + AS4 + AS5 + AS6 + AS7 + AS8 + AS9 + AS10"
m2 <- "SC =~ SC1 + SC2 + SC3 + SC4 + SC5 + SC6 + SC7 + SC8 + SC9 + SC10"
m3 <- "AD =~ AD1 + AD2 + AD3 + AD4 + AD5 + AD6 + AD7 + AD8 + AD9 + AD10"
m4 <- "DO =~ DO1 + DO2 + DO3 + DO4 + DO5 + DO6 + DO7 + DO8 + DO9 + DO10"

getLat <- function(x){
  fit <- sem(x, df, estimator="wlsmv", ordered=names(df)[1:40], group = "gender", group.equal=c("thresholds", "loadings", "intercepts"), std.lv=F, auto.fix.first=T)
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

lv <- list(getLat(m1), getLat(m2), getLat(m3), getLat(m4))
lat <- unlist(lapply(lv, `[[`, 1))

#obs diff
for (i in 1:40) {
  df[[i]] <- as.numeric(paste(df[[i]]))
}
df$AS_sum <- scale(rowSums(df[c("AS1", "AS2", "AS3", "AS4", "AS5", "AS6", "AS7", "AS8", "AS9", "AS10")]))
df$SC_sum <- scale(rowSums(df[c("SC1", "SC2", "SC3", "SC4", "SC5", "SC6", "SC7", "SC8", "SC9", "SC10")]))
df$AD_sum <- scale(rowSums(df[c("AD1", "AD2", "AD3", "AD4", "AD5", "AD6", "AD7", "AD8", "AD9", "AD10")]))
df$DO_sum <- scale(rowSums(df[c("DO1", "DO2", "DO3", "DO4", "DO5", "DO6", "DO7", "DO8", "DO9", "DO10")]))

get_d <- function(x){
  v1 <- x[df$gender==1] 
  v2 <- x[df$gender==2]
  mean_diff <- mean(v1)-mean(v2)
  pool_sd <- sqrt((var(v1)+var(v2))/2)
  d <- mean_diff/pool_sd
  return(d)
}

obs <- c(get_d(df$AS_sum),
         get_d(df$SC_sum),
         get_d(df$AD_sum),
         get_d(df$DO_sum))


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
