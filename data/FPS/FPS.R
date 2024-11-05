#observed score vs latent score, gender diffs
library(lavaan)

df <- read.csv("data.csv", sep="\t")
df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[,1:60]<1)==0,]

for (i in 1:60) {
  df[[i]] <- ordered(df[[i]], levels=sort(unique(df[[i]])))
}

#recode items
m1 <- "FPS_Conservative =~ Q1 + Q4 + Q13 + Q17 + Q23 + Q36 + Q38 + Q47 + Q53 + Q59"
m2 <- "FPS_Liberal =~ Q5 + Q6 + Q7 + Q22 + Q24 + Q27 + Q33 + Q42 + Q52 + Q60"
m3 <- "FPS_Radical =~ Q2 + Q15 + Q16 + Q18 + Q19 + Q29 + Q34 + Q46 + Q48 + Q55"
m4 <- "FPS_Socialist =~ Q10 + Q20 + Q25 + Q31 + Q39 + Q41 + Q45 + Q54 + Q56 + Q58"
m5 <- "FPS_Cultural =~ Q9 + Q11 + Q14 + Q28 + Q30 + Q32 + Q35 + Q37 + Q44 + Q50"
m6 <- "FPS_Color =~ Q3 + Q8 + Q12 + Q21 + Q26 + Q40 + Q43 + Q49 + Q51 + Q57"

#lat diff
getLat <- function(x){
  fit <- sem(x, df, estimator="wlsmv", ordered=names(df)[1:60], group = "gender", group.equal=c("thresholds", "loadings", "intercepts"), std.lv=F, auto.fix.first=T)
  omega <- 0
  omega <- mean(semTools::compRelSEM(fit, ord.scale = T))
  if (is.na(omega)) {
    l1<-lavInspect(fit, what="std")$`1`$lambda; l2<-lavInspect(fit, what="std")$`2`$lambda
    omega <- mean(sum(abs(l1))^2 / (sum(abs(l1))^2 + sum(1-abs(l1)^2)), sum(abs(l2))^2 / (sum(abs(l2))^2 + sum(1-abs(l2^2))))
  }
  fitM <- fitmeasures(fit, fit.measures = c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  return(list(lavInspect(fit, what="std")$`2`$alpha, 
              omega,
              fitM))
}

lv <- list(getLat(m1), getLat(m2), getLat(m3), getLat(m4), getLat(m5), getLat(m6))
lat <- -1*unlist(lapply(lv, `[[`, 1))


#obs diff
for (i in 1:60) {
  df[[i]] <- as.numeric(paste(df[[i]]))
}
df$FPS_Conservative <- scale(rowSums(df[,c(1,4,13,17,23,36,38,47,53,59)]))
df$FPS_Liberal <- scale(rowSums(df[,c(5,6,7,22,24,27,33,42,52,60)]))
df$FPS_Radical <- scale(rowSums(df[,c(2,15,16,18,19,29,34,46,48,55)]))
df$FPS_Socialist <- scale(rowSums(df[,c(10,20,25,31,39,41,45,54,56,58)]))
df$FPS_Cultural <- scale(rowSums(df[,c(9,11,14,28,30,32,35,37,44,50)]))
df$FPS_Color <- scale(rowSums(df[,c(3,8,12,21,26,40,43,49,51,57)]))


get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- c(get_d(df$FPS_Conservative),
         get_d(df$FPS_Liberal),
         get_d(df$FPS_Radical),
         get_d(df$FPS_Socialist),
         get_d(df$FPS_Cultural),
         get_d(df$FPS_Color))

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