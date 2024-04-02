library(Rcpp)
library(pracma)
library(GPArotation)
library(MASS)
load("~/Desktop/code_new/data_result.RData")
load("~/Desktop/data_analysis_new.RData")
K <- 3
grid <- 51

true_rate <- array(rep(0,grid*n_individual*n_item), dim=c(grid,n_individual,n_item))
for(k in 1:grid){
  true_rate[k, , ] <- Rate_left[k, , ] %*% Rate_right
}

Rate_left_sum <- matrix(0,nrow = n_individual, ncol = K)
for(k in 1:grid){
  Rate_left_sum <- Rate_left_sum + Rate_left[k,,]/51
}
Rate_mean <- colMeans(Rate_left_sum)
Rate_left_sum <- Rate_left_sum - t(matrix(rep(Rate_mean,nrow(Rate_left_sum)),nrow = 3))

for(k in 1:grid){
  Rate_left[k,,] <- Rate_left[k,,] - t(matrix(rep(Rate_mean,nrow(Rate_left_sum)),nrow = 3))
}

d <- eigen(t(Rate_left_sum)%*%Rate_left_sum)
R <- d$vectors %*% diag((d$values)^(1/2)) %*% t(d$vectors)

A <- t(R %*% Rate_right)
B <- varimax(A, normalize = T)

rate_right_rotated <- B$loadings
rank <- vector("list", K)
for(i in 1:K){
  Order <- order(abs(rate_right_rotated[,i]), decreasing=T)
  rank[[i]] <- data.frame(value = rate_right_rotated[Order, i],
                          product = product_detail$COMMODITY_DESC[Order],
                          sub = product_detail$SUB_COMMODITY_DESC[Order])
  rank[[i]] <- rank[[i]][(rank[[i]]$product!=" ")&(rank[[i]]$sub!="NO SUBCOMMODITY DESCRIPTION"),]
}

rank_cross <- vector("list", K)
for(i in 1:K){
  cross <- rep(0, n_item)
  cross <- rate_right_rotated[,i]/sqrt(rate_right_rotated[,1]^2+rate_right_rotated[,2]^2+rate_right_rotated[,3]^2)
  Order <- order(cross, decreasing = T)
  rank_cross[[i]] <- data.frame(value = cross[Order],
                          product = product_detail$COMMODITY_DESC[Order],
                          sub = product_detail$SUB_COMMODITY_DESC[Order])
  rank_cross[[i]] <- rank_cross[[i]][(rank_cross[[i]]$product!=" ")&(rank_cross[[i]]$sub!="NO SUBCOMMODITY DESCRIPTION"),]
}

a_1 <- rank_cross[[1]][1:60,]
a_2 <- rank_cross[[2]][1:60,]
a_3 <- rank_cross[[3]][1:60,]
prod_1 <- sort(table(rank_cross[[1]]$sub[1:50]), decreasing = TRUE)
prod_2 <- sort(table(rank_cross[[2]]$sub[1:50]), decreasing = TRUE)
prod_3 <- sort(table(rank_cross[[3]]$sub[1:50]), decreasing = TRUE)
prod_1[prod_1 > 1]
prod_2[prod_2 > 1]
prod_3[prod_3 > 1]

rate_left_rotated <- array(rep(0,grid * n_individual * K), dim=c(grid, n_individual, K))
for(k in 1:grid){
  rate_left_rotated[k,,] <- Rate_left[k,,] %*% solve(R) %*% B$rotmat %*% diag(c(1,1,1))
}

flag <- match(individual_info$household_key,individual_id)
individual_info <- individual_info[is.na(flag) == 0,]
flag <- flag[is.na(flag) == 0]
individual_info$AGE_DESC <- as.numeric(factor(individual_info$AGE_DESC))
individual_info$INCOME_DESC <- as.numeric(factor(individual_info$INCOME_DESC,
                                                 levels = c("Under 15K", "15-24K", "25-34K", "35-49K",
                                                            "50-74K", "75-99K", "100-124K", "125-149K", 
                                                            "150-174K", "175-199K", "200-249K", "250K+"), 
                                                 labels = 1:12))
individual_info$HOUSEHOLD_SIZE_DESC <- as.numeric(factor(individual_info$HOUSEHOLD_SIZE_DESC))
individual_info$KID_CATEGORY_DESC <- as.numeric(factor(individual_info$KID_CATEGORY_DESC,
                                                 levels = c("None/Unknown", "1", "2", "3+"), 
                                                 labels = 1:4))
individual_info <- individual_info[,c("AGE_DESC", "INCOME_DESC", "HOUSEHOLD_SIZE_DESC", "KID_CATEGORY_DESC")]
individual_info$KID_CATEGORY_DESC[individual_info$KID_CATEGORY_DESC==1] <- 0
individual_info$KID_CATEGORY_DESC[individual_info$KID_CATEGORY_DESC>1] <- 1
individual_info$AGE_DESC[individual_info$AGE_DESC <= 4] <- 0
individual_info$AGE_DESC[individual_info$AGE_DESC > 1] <- 1
individual_info$INCOME_DESC[individual_info$INCOME_DESC <= 3] <- 0
individual_info$INCOME_DESC[individual_info$INCOME_DESC <= 5&individual_info$INCOME_DESC > 0] <- 1
individual_info$INCOME_DESC[individual_info$INCOME_DESC > 1] <- 2
individual_info <- individual_info[,-3]

rate_left_rotated <- rate_left_rotated[,,c(3,2,1)]
individual_info$factor_1 <- rep(0,nrow(individual_info))
individual_info$factor_2 <- rep(0,nrow(individual_info))
individual_info$factor_3 <- rep(0,nrow(individual_info))
for(i in 1:nrow(individual_info)){
  for(j in 1:3){
    individual_info[i,(j+3)] <- mean(rate_left_rotated[,flag[i],j])
  }
}

cov_matrix <- matrix(0, 6, 6)
cov_test <- matrix(0, 6, 6)
for (i in 1:6) {
  for (j in 1:6) {
    if(i>3 & j>3){
      cov_matrix[i, j] <- cor(individual_info[, i], individual_info[, j])
      test <- cor.test(individual_info[, i], individual_info[, j])
      cov_test[i, j] <- test$p.value
    }else{
      cov_matrix[i, j] <- cor(individual_info[, i], individual_info[, j])
      test <- cor.test(individual_info[, i], individual_info[, j])
      cov_test[i, j] <- test$p.value
    }
  }
}

TV <- function(x){
  l <- length(x)
  sum <- 0
  for(i in 2:l){
    sum <- sum + abs(x[i] - x[i-1])/(l-1)
  }
  return(sum)
}

var_total <- matrix(rep(0,n_item*n_individual),nrow=n_item)
for(j in 1:n_item){
  record <- rep(0,n_individual)
  for(i in 1:n_individual){
    record[i] <- TV(true_rate[,i,j])
    if(record[i]==0){
      print(i)
    }
  }
  var_total[j,] <- record
}

transaction$type <- product_detail$COMMODITY_DESC[transaction$PRODUCT_ID]
product_number <- 18
type <- names(sort(table(transaction$type), decreasing = TRUE)[1:product_number])

fluc <- data.frame(type = type, quant_1 = rep(0, length(type)), quant_2 = rep(0, length(type)), quant_3 = rep(0, length(type)))
for(i in 1:length(type)){
  fluc$quant_1[i] <- quantile(var_total[product_detail$COMMODITY_DESC == type[i],], probs = 0.25)
  fluc$quant_2[i] <- quantile(var_total[product_detail$COMMODITY_DESC == type[i],], probs = 0.5)
  fluc$quant_3[i] <- quantile(var_total[product_detail$COMMODITY_DESC == type[i],], probs = 0.75)
}

fluc_1 <- fluc[order(fluc$quant_1, decreasing = T),]
fluc_2 <- fluc[order(fluc$quant_2, decreasing = T),]
fluc_3 <- fluc[order(fluc$quant_3, decreasing = T),]

par(mfrow = c(1, 3))
par(mar=c(20, 5, 2, 2))
barplot(fluc_1$quant_1,names.arg = fluc_1$type, cex.axis = 2, cex.main = 2, cex.names = 1.5, las = 2, main = "First Quatile")
barplot(fluc_2$quant_2,names.arg = fluc_2$type, cex.axis = 2, cex.main = 2, cex.names = 1.5, las = 2, main = "Median")
barplot(fluc_3$quant_3,names.arg = fluc_3$type, cex.axis = 2, cex.main = 2, cex.names = 1.5, las = 2, main = "Third Quatile")

lm_1 = lm(factor_1 ~ AGE_DESC+factor(INCOME_DESC)+KID_CATEGORY_DESC1+
         +factor(INCOME_DESC*KID_CATEGORY_DESC), data=individual_info)
summary(lm_1)

lm_2 = lm(factor_2 ~ AGE_DESC+factor(INCOME_DESC)+KID_CATEGORY_DESC
         +factor(INCOME_DESC*KID_CATEGORY_DESC), data=individual_info)
summary(lm_2)

lm_3 = lm(factor_3 ~ AGE_DESC+factor(INCOME_DESC)+KID_CATEGORY_DESC
         +factor(INCOME_DESC*KID_CATEGORY_DESC), data=individual_info)
summary(lm_3)
