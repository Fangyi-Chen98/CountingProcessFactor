product <- read.csv("~/Desktop/HDCP/data/product.csv")
individual_info <- read.csv("~/Desktop/HDCP/data/hh_demographic.csv")
transaction <- read.csv("~/Desktop/HDCP/data/transaction_data.csv")
transaction <- transaction[,c("household_key", "DAY", "PRODUCT_ID", "TRANS_TIME")]
transaction$hour <- floor(transaction$TRANS_TIME/100)
transaction$minute <- transaction$TRANS_TIME - transaction$hour * 100
transaction$time <- (transaction$DAY - 1) * 1440 + transaction$hour * 60 + transaction$minute
transaction <- transaction[, c("household_key", "PRODUCT_ID", "time")]
transaction$time <- transaction$time/(711*24*60)
transaction <- transaction[transaction$time>0.15,]
transaction$time <- (transaction$time - 0.15) / (1 - 0.15)

transaction <- transaction[,c("household_key", "PRODUCT_ID")]
label <- match(transaction$PRODUCT_ID, product$PRODUCT_ID)
transaction$COMMODITY_DESC <- product$COMMODITY_DESC[label]
transaction$SUB_COMMODITY_DESC <- product$SUB_COMMODITY_DESC[label]

product_number <- 2000
most <- as.numeric(names(sort(table(transaction$PRODUCT_ID), decreasing = TRUE)[1:product_number]))
transaction <- transaction[transaction$PRODUCT_ID %in% most == 1,]
rm(most)

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

rm(number_individual)
rm(number_product)
rm(match_individual)
rm(match_product)
transaction <- transaction[,-1]

label <- match(1:product_number, transaction$PRODUCT_ID)
product_detail <- transaction[label,]
rm(flag)
rm(n_individual)
rm(label)
rm(transaction)
rm(product)
rm(n_item)
save.image("~/Desktop/code_new/data_result.RData")