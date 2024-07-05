#observed score vs latent score, gender diffs
library(lavaan)

df <- read.csv("data.csv",sep = "\t")
df <- df[c(grep("A",names(df)),grep("gender",names(df)))]

df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[,1:56]==0)==0,]


#recode items, none
model <- "FTI_CURIOUS =~ Q1A + Q2A + Q3A + Q4A + Q5A + Q6A + Q7A + Q8A + Q9A + Q10A + Q11A + Q12A + Q13A + Q14A
          FTI_CAUTIOUS =~ Q15A + Q16A + Q17A + Q18A + Q19A + Q20A + Q21A + Q22A + Q23A + Q24A + Q25A + Q26A + Q27A + Q28A
          FTI_ANALYTICAL =~ Q29A + Q30A + Q31A + Q32A + Q33A + Q34A + Q35A + Q36A + Q37A + Q38A + Q39A + Q40A + Q41A + Q42A
          FTI_PROSOCIAL =~ Q43A + Q44A + Q45A + Q46A + Q47A + Q48A + Q49A + Q50A + Q51A + Q52A + Q53A + Q54A + Q55A + Q56A"

fit <- sem(model, df, estimator="mlr", group = "gender", group.equal=c("loadings", "intercepts"), std.lv=T, auto.fix.first=F)
lat <- (lavInspect(fit, what="std")$`1`$alpha)


#obs diff
df$CURIOUS_sum <- scale(rowSums(df[1:14]))
df$CAUTIOUS_sum <- scale(rowSums(df[15:28]))
df$ANALYTICAL_sum <- scale(rowSums(df[29:42]))
df$PROSOCIAL_sum <- scale(rowSums(df[43:56]))

get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- c(get_d(df$CURIOUS_sum),
         get_d(df$CAUTIOUS_sum),
         get_d(df$ANALYTICAL_sum),
         get_d(df$PROSOCIAL_sum))

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