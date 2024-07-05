#observed score vs latent score, gender diffs

library(lavaan)

df <- read.csv("data.csv", sep="\t")

temp <- df[1:126]
df <- cbind(temp[grep("A", names(temp))],
            df[c("gender")])
rm(temp)

df <- df[df$gender %in% c(1,2),]
df <- df[rowSums(df[,1:42]<1)==0,]
table(df$Q1A)

#recode items <- no inverted items
model <- "DASS_DEP =~ Q3A + Q5A + Q10A + Q13A + Q16A + Q17A + Q21A + Q24A + Q26A + Q31A + Q34A + Q37A + Q38A + Q42A
          DASS_ANX =~ Q2A + Q4A + Q7A + Q9A + Q15A + Q19A + Q20A + Q23A + Q25A + Q28A + Q30A + Q36A + Q40A + Q41A
          DASS_STR =~ Q1A + Q6A + Q8A + Q11A + Q12A + Q14A + Q18A + Q22A + Q27A + Q29A + Q32A + Q33A + Q35A + Q39A"

fit <- sem(model, df, estimator="mlr", group = "gender", group.equal=c("loadings", "intercepts"), std.lv=T, auto.fix.first=F)
lat <- (lavInspect(fit, what="std")$`1`$alpha)


#obs diff
df$DEP_sum <- scale(rowSums(df[,c(3,5,10,13,16,17,21,24,26,31,34,37,38,42)]))
df$ANX_sum <- scale(rowSums(df[,c(2,4,7,9,15,19,20,23,25,28,30,36,40,41)]))
df$STR_sum <- scale(rowSums(df[,c(1,6,8,11,12,14,18,22,27,29,32,33,35,39)]))

get_d <- function(x){
    v1 <- x[df$gender==1] 
    v2 <- x[df$gender==2]
    mean_diff <- mean(v1)-mean(v2)
    pool_sd <- sqrt((var(v1)+var(v2))/2)
    d <- mean_diff/pool_sd
    return(d)
}

obs <- c(get_d(df$DEP_sum),
         get_d(df$ANX_sum),
         get_d(df$STR_sum))

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
