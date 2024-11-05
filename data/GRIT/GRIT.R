#observed score vs latent score, gender diffs
library(lavaan)

df <- read.csv("data.csv", sep="\t")
df <- df[df$gender %in% c(1,2),]

temp <- df[3:14]
df <- cbind(temp, df["gender"])

table(df$GS1)
df <- df[rowSums(df[,1:10]<1)==0,]

#recode items
m1 <- "GRIT =~ GS1 + GS2 + GS3 + GS4 + GS5 + GS6 + GS7 + GS8 + GS9 + GS10 + GS11 + GS12"

rec_vec <- 1:12
rec_vec <- rec_vec[lavInspect(sem(m1, df), what="std")$lambda>0]
for (i in rec_vec) {
    df[,i] <- (df[,i]-6)*(-1)
}

for (i in 1:12) {
  df[[i]] <- ordered(df[[i]], levels=sort(unique(df[[i]])))
}

#lat diff
getLat <- function(x){
  fit <- sem(x, df, estimator="wlsmv", ordered=names(df)[1:12], group = "gender", group.equal=c("thresholds", "loadings", "intercepts"), std.lv=F, auto.fix.first=T)
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

lv <- list(getLat(m1))
lat <- unlist(lapply(lv, `[[`, 1))

#obs diff
for (i in 1:12) {
  df[[i]] <- as.numeric(paste(df[[i]]))
}
df$GS_sum <- scale(rowSums(df[,1:12]))

get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- get_d(df$GS_sum)

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
