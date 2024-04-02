.libPaths("/moto/stats/users/fc2630/rpackages")
library(Rcpp)
K <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(K)
load("initialization_new.RData")
Rcpp::sourceCpp('data_analysis.cpp')
grid_num <- 50
n_individual <- 1978
n_item <- 2000
max_iter <- 600
max_grad <- 40
step_size_initial <- 1
h <- 0.1*(min(c(n_item,n_individual))/(log(min(c(n_item,n_individual)))*log(min(c(n_item,n_individual)))))^(-1/5)
t <- seq(h,1-h,(1-2*h)/grid_num)
M <- 4
threshold <- 0.025

rate_initial <- vector("list", length = length(t)+1)
for(l in 1:length(t)){
  rate <- Rate_left[[K]]
  rate_initial[[l]] <- rate[l,,]
}
rate_initial[[length(t)+1]] <- Rate_right[[K]]
result <- data_analysis(as.matrix(transaction), nrow(transaction), rate_initial, n_individual, n_item, K, max_grad, max_iter, step_size_initial, M, threshold, grid_num + 1, min(t), max(t), h)
if(K > 1){
  rate_log_left <- array(rep((grid_num + 1)*n_individual*K),dim=c((grid_num + 1), n_individual, K))
  rate_log_right <- result[[grid_num + 2]]
  for(l in 1:length(t)){
    rate_log_left[l,,] <- result[[l]]
  }
  Rate_left <- rate_log_left
  Rate_right <- rate_log_right
}else{
  rate_log_left <- array(rep((grid_num + 1)*n_individual*2),dim=c((grid_num + 1), n_individual, 2))
  rate_log_right <- result[[grid_num + 2]]
  for(l in 1:length(t)){
    rate_log_left[l,,] <- result[[l]]
  }
  Rate_left <- rate_log_left
  Rate_right <- rate_log_right
}
Log_lik_int_ker <- result[[grid_num + 3]]
rm(result)
print(Log_lik_int_ker)
save.image(paste0("data_analysis_new_",as.character(K),".RData"))