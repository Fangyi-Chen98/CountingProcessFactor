library(Rcpp)
setwd("~/Desktop/HDCP/")
product <- read.csv("data/product.csv")
transaction <- read.csv("data/transaction_data.csv")
transaction <- transaction[,c("household_key", "DAY", "PRODUCT_ID", "TRANS_TIME")]
transaction$hour <- floor(transaction$TRANS_TIME/100)
transaction$minute <- transaction$TRANS_TIME - transaction$hour * 100
transaction$time <- (transaction$DAY - 1) * 1440 + transaction$hour * 60 + transaction$minute
transaction <- transaction[, c("household_key", "PRODUCT_ID", "time")]
transaction$time <- transaction$time/(711*24*60)
transaction <- transaction[transaction$time>0.15,]
transaction$time <- (transaction$time - 0.15) / (1 - 0.15)

product_number <- 2000
most <- as.numeric(names(sort(table(transaction$PRODUCT_ID), decreasing = TRUE)[1:product_number]))
transaction <- transaction[transaction$PRODUCT_ID %in% most == 1,]
rm(most)
rm(product)
id <- unique(transaction$household_key)
flag <- rep(0, length(id))
for(i in 1:length(id)){
  if(sum(transaction$household_key == id[i]) >= 100){
    flag[i] <- 1
  }
}
id <- id[flag == 1]

Flag <- rep(0, nrow(transaction))
for(i in 1:nrow(transaction)){
  Flag[i] <- transaction$household_key[i] %in% id
}
transaction <- transaction[Flag == 1, ]

product_id <- sort(unique(transaction$PRODUCT_ID))
match_product <- function(id){
  return(match(id,product_id))
}
number_product <- sapply(transaction$PRODUCT_ID, match_product)
individual_id <- sort(unique(transaction$household_key))
match_individual <- function(id){
  return(match(id,individual_id))
}
number_individual <- sapply(transaction$household_key, match_individual)
transaction$household_key <- number_individual
transaction$PRODUCT_ID <- number_product
n_individual <- length(individual_id)
n_item <- length(product_id)
rm(individual_id)
rm(number_individual)
rm(number_product)
rm(product_id)
rm(match_individual)
rm(match_product)

Rcpp::sourceCpp('code/C/count_event.cpp')
K_max <- 5
grid_num <- 50
count <- matrix(rep(0,n_individual*n_item),nrow=n_individual)
count_event(count,as.matrix(transaction),nrow(transaction),0,1)
rate_int <- log(count+0.01)
rate_int_svd <- svd(rate_int)
Rate_left <- vector("list",K_max)
Rate_right <- vector("list",K_max)
for(k in 1:K_max){
  if(k == 1){
    rate <- array(rep((grid_num + 1)*n_individual*2),dim=c((grid_num + 1), n_individual, 2))
    d <- rate_int_svd$d
    rate_left <- rate_int_svd$u[,1:2]%*% sqrt(diag(c(d[1],0)))
    Rate_right[[k]] <- sqrt(diag(c(d[1],0)))%*% t(rate_int_svd$v[,1:2])
    for(l in 1:(grid_num + 1)){
      rate[l,,] <- rate_left
    }
    Rate_left[[k]] <- rate
  }else{
    rate <- array(rep((grid_num + 1)*n_individual*k),dim=c((grid_num + 1), n_individual, k))
    d <- rate_int_svd$d
    rate_left <- rate_int_svd$u[,1:k]%*% sqrt(diag(d[1:k]))
    Rate_right[[k]] <- sqrt(diag(d[1:k]))%*% t(rate_int_svd$v[,1:k])
    for(l in 1:(grid_num + 1)){
      rate[l,,] <- rate_left
    }
    Rate_left[[k]] <- rate
  }
}

setwd("~/Desktop/")
save.image("code_new/initialization_new.RData")
