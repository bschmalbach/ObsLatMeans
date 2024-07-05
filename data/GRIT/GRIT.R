#observed score vs latent score, gender diffs
library(lavaan)

df <- read.csv("data.csv", sep="\t")
df <- df[df$gender %in% c(1,2),]

temp <- df[3:14]
df <- cbind(temp, df["gender"])

table(df$GS1)
df <- df[rowSums(df[,1:10]<1)==0,]

#recode items
model <- "GRIT =~ GS1 + GS2 + GS3 + GS4 + GS5 + GS6 + GS7 + GS8 + GS9 + GS10 + GS11 + GS12"

rec_vec <- 1:12
rec_vec <- rec_vec[lavInspect(sem(model, df), what="std")$lambda>0]
for (i in rec_vec) {
    df[,i] <- (df[,i]-6)*(-1)
}


#lat diff
fit <- sem(model, df, estimator="mlr", group = "gender", group.equal=c("loadings", "intercepts"), std.lv=T, auto.fix.first=F)
lat <- (lavInspect(fit, what="std")$`1`$alpha)


#obs diff
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

res$z <- res$diff/(res$pool_sd/sqrt(sum(n_g)))
res$p <- 2*(1-pnorm(abs(res$z)))
res$dd <- res$diff/(res$pool_sd)

res$UL_CI <- res$diff + 1.96*res$pool_sd
res$LL_CI <- res$diff - 1.96*res$pool_sd
res$n <- sum(n_g)
res$omega <-  rowMeans(matrix(unlist(semTools::reliability(fit, what="omega")),nrow = nrow(res), byrow = F))
