#observed score vs latent score, gender diffs
library(lavaan)

df <- read.csv("data.csv",sep = "\t")
df <- df[df$gender %in% c(1,2),]
df <- df[grepl("A", names(df)) | grepl("gender", names(df))]
df <- df[rowSums(df[1:19]==0)==0,]

#recode items
rec_vec <- c(9, 13, 15)
for (i in rec_vec) {
  df[,i] <- (6-df[,i])
}

for (i in 1:19) {
  df[[i]] <- ordered(df[[i]], levels=sort(unique(df[[i]])))
}

m1 <- "PWE =~ Q1A + Q2A + Q3A + Q4A + Q5A + Q6A + Q7A + Q8A + Q9A + Q10A + Q11A + Q12A + Q13A + Q14A + Q15A + Q16A + Q17A + Q18A + Q19A"

getLat <- function(x){
  fit <- sem(x, df, estimator="wlsmv", ordered=names(df)[1:19], group = "gender", group.equal=c("thresholds", "loadings", "intercepts"), std.lv=F, auto.fix.first=T)
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

lv <- list(getLat(m1))
lat <- unlist(lapply(lv, `[[`, 1))


#obs diff
for (i in 1:19) {
  df[[i]] <- as.numeric(paste(df[[i]]))
}
df$PWE_sum <- scale(rowSums(df[1:19]))

get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- get_d(df$PWE_sum)

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