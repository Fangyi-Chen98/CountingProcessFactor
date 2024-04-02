.libPaths("/moto/stats/users/fc2630/rpackages")
library(Rcpp)
library(pracma)
id <- as.character(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(as.numeric(id))
print(as.numeric(id))
myFunction <- cppFunction("
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
//[[Rcpp::export]]
int generate_data_constant(NumericMatrix event, int n_individual, int n_item, int K, int max_event, double mu1, double mu2, double mu3, double mu4, double mu5, double mu6,NumericMatrix Rate, NumericMatrix Rate_left, NumericMatrix Rate_right){
    GetRNGstate();
    int i,j,m,event_count=0;
    double time_left, sum, current_time,T;
    double X_left[n_individual][K],X_right[K][n_item];
    T = 1.00;
    for(i=0;i<n_individual;i++){
        for(j=0;j<K;j++){
            double E;
            E = R::runif(mu1, mu2);
            X_left[i][j] = E;
            Rate_left(i,j) = E;
        }
    }
    for(i=0;i<K-1;i++){
        for(j=0;j<n_item;j++){
            double E;
            E = R::runif(mu3,mu4);
            X_right[i][j] = E;
            Rate_right(i,j) = E;
        }
    }
    for(j=0;j<n_item;j++){
        double E;
        E = R::runif(mu5,mu6);
        X_right[K-1][j] = E;
        Rate_right((K-1),j) = E;
    }
    for(i=0;i<n_individual;i++){
        for(j=0;j<n_item;j++){
            sum = 0.0;
            for(m=0; m<K; m++){
                sum = sum + X_left[i][m] * X_right[m][j];
            }
            Rate(i,j) = sum;
        }
    }
    double rate;
    for(i=0;i<n_individual;i++){
        for(j=0;j<n_item;j++){
            rate = exp(Rate(i,j));
            time_left = rate * T;
            current_time = 0.0;
            while(time_left>0){
                double U = R::rexp(1.0);
                if(time_left > U){
                    current_time = current_time + U;
                    event(event_count,0) = i+1;
                    event(event_count,1) = j+1;
                    event(event_count,2) = current_time / rate;
                    event_count++;
                    time_left = time_left - U;
                    if(event_count == max_event){
                        return -1;
                    }
                }else{
                    break;
                }
            }
        }
    }
    return event_count;
}
")
myFunction <- cppFunction("
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
//[[Rcpp::export]]
int generate_data_linear(NumericMatrix event, int n_individual, int n_item, int K, int max_event, double mu1, double mu2, double mu3, double mu4, double mu5, double mu6, NumericMatrix Rate, NumericMatrix Rate_left, NumericMatrix Rate_right, NumericMatrix coef){
    GetRNGstate();
    int i,j,m,event_count=0;
    double time_left, sum, current_time,T;
    double X_left[n_individual][K],X_right[K][n_item];
    T = 1.00;
    for(i=0;i<n_individual;i++){
        for(j=0;j<K;j++){
            double E;
            E = R::runif(mu1,mu2);
            X_left[i][j] = E-0.5*coef(i,j);
            Rate_left(i,j) = E-0.5*coef(i,j);
        }
    }
    for(i=0;i<K-1;i++){
        for(j=0;j<n_item;j++){
            double E;
            E = R::runif(mu3,mu4);
            X_right[i][j] = E;
            Rate_right(i,j) = E;
        }
    }
    for(j=0;j<n_item;j++){
        double E;
        E = R::runif(mu5,mu6);
        X_right[K-1][j] = E;
        Rate_right((K-1),j) = E;
    }
    for(i=0;i<n_individual;i++){
        for(j=0;j<n_item;j++){
            sum = 0.0;
            for(m=0; m<K; m++){
                sum = sum + X_left[i][m] * X_right[m][j];
            }
            Rate(i,j) = sum;
        }
    }
    double rate,a,b,proposed_time;
    for(i=0;i<n_individual;i++){
        for(j=0;j<n_item;j++){
            a = Rate(i,j);
            rate = exp(a);
            b = 0;
            for(m=0;m<K;m++){
                b += coef(i,m) * X_right[m][j];
            }
            time_left = (rate*exp(b*T)-rate)/b;
            current_time = 0.0;
            while(time_left>0){
                double U = R::rexp(1.0);
                if(time_left > U){
                    current_time = current_time + U;
                    event(event_count,0) = i+1;
                    event(event_count,1) = j+1;
                    event(event_count,2) = (log(b*current_time+exp(a))-a)/b;
                    event_count++;
                    time_left = time_left - U;
                    if(event_count == max_event){
                        return -1;
                    }
                }else{
                    break;
                }
            }
        }
    }
    return event_count;
}
")
myFunction <- cppFunction("
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
//[[Rcpp::export]]
int generate_data_periodic(NumericMatrix event, int n_individual, int n_item, int K, int max_event, double mu1, double mu2, double mu3, double mu4, double mu5, double mu6, NumericMatrix Rate, NumericMatrix Rate_left,NumericMatrix Rate_right, NumericMatrix coef_scale, NumericMatrix coef_phase, double period){
    GetRNGstate();
    int i,j,m,event_count=0;
    double time_left, sum, current_time,T;
    double X_left[n_individual][K],X_right[K][n_item];
    T = 1.00;
    for(i=0;i<n_individual;i++){
        for(j=0;j<K;j++){
            double E;
            E = R::runif(mu1,mu2);
            X_left[i][j] = E;
            Rate_left(i,j) = E;
        }
    }
    for(i=0;i<K-1;i++){
        for(j=0;j<n_item;j++){
            double E;
            E = R::runif(mu3,mu4);
            X_right[i][j] = E;
            Rate_right(i,j) = E;
        }
    }
    for(j=0;j<n_item;j++){
        double E;
        E = R::runif(mu5,mu6);
        X_right[K-1][j] = E;
        Rate_right((K-1),j) = E;
    }
    for(i=0;i<n_individual;i++){
        for(j=0;j<n_item;j++){
            sum = 0.0;
            for(m=0; m<K; m++){
                sum = sum + X_left[i][m] * X_right[m][j];
            }
            Rate(i,j) = sum;
        }
    }
    double rate_max, rate_current;
    for(i=0;i<n_individual;i++){
        for(j=0;j<n_item;j++){
            rate_max = Rate(i,j);
            for(m=0; m<K; m++){
                rate_max += fabs(coef_scale(i,m)*X_right[m][j]);
            }
            time_left = T;
            current_time = 0.0;
            while(time_left>0){
                double U = R::rexp(exp(-rate_max));
                if(time_left > U){
                    current_time = current_time + U;
                    time_left = time_left - U;
                    rate_current = Rate(i,j);
                    for(m=0; m<K; m++){
                        rate_current += (coef_scale(i,m)*X_right[m][j])*sin(2*M_PI*(current_time+coef_phase(i,m))/period);
                    }
                    double U_1 = R::runif(0,1);
                    if(U_1 < exp(rate_current-rate_max)){
                        event(event_count,0) = i+1;
                        event(event_count,1) = j+1;
                        event(event_count,2) = current_time;
                        event_count++;
                    }
                    if(event_count == max_event){
                        return -1;
                    }
                }else{
                    break;
                }
            }
        }
    }
    return event_count;
}
")
myFunction <- cppFunction("
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
//[[Rcpp::export]]
List poisson_constant(NumericMatrix event, int nrow, NumericMatrix T_mat, NumericMatrix rate_left_true, NumericMatrix rate_right_true, int row, int col, int K, int K_true, int max_grad, int max_iter, double step_size_initial, double M, double threshold, int grid_num, double t_min, double t_max) {
	  int i, j, k, count_iter, count_grad, row_event, col_event;
    double grad, hess, sum, f_new = 0.0, f_old = 0.0, step_size, grid_size;
    double product[row][col];
    double weight[row][col];
    double B_old[K][col], B_new[K][col], A_new[row][K], A_old[row][K];
    grid_size = (t_max - t_min)/(grid_num - 1);
    double A[row][K], B[K][col];
    if(K <= K_true){
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_new[i][j] = rate_left_true(i,j);
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = rate_right_true(i,j);
            }
        }
    }else{
        for(i=0; i<row; i++){
            for(j=0; j<K_true; j++){
                A[i][j] = rate_left_true(i,j);
            }
            for(j=K_true; j<K; j++){
                A[i][j] = 0.0;
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_new[i][j] = 0.0;
                for(k=0; k<K; k++){
                    A_new[i][j] += A[i][k] * T_mat(k,j);
                }
            }
        }
        for(i=0; i<K_true; i++){
            for(j=0; j<col; j++){
                B[i][j] = rate_right_true(i,j);
            }
        }
        for(i=K_true; i<K; i++){
            for(j=0; j<col; j++){
                B[i][j] = 0.0;
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = 0.0;
                for(k=0; k<K; k++){
                    B_new[i][j] += T_mat(k,i) * B[k][j];
                }
            }
        }
    }
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            weight[i][j] = 0.0;
        }
    }
    for(i=0;i<nrow;i++){
        row_event = event(i,0) - 1;
        col_event = event(i,1) - 1;
        weight[row_event][col_event] ++;
    }
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            product[i][j] = 0.0;
            for(k=0; k<K; k++){
                product[i][j] = product[i][j] + A_new[i][k] * B_new[k][j];
            }
            f_new = f_new - (weight[i][j] * product[i][j] - exp(product[i][j]));
        }
    }
    count_iter = 0;
    while(count_iter == 0|| fabs(f_new - f_old) > threshold){
        f_old = f_new;
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_old[i][j] = B_new[i][j];
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_old[i][j] = A_new[i][j];
            }
        }
        count_grad = 0;
        step_size = step_size_initial;
        while(f_new >= f_old){
            step_size = step_size * 0.5;
            if(count_grad > 0){
                for(i=0; i<K; i++){
                    for(j=0; j<col; j++){
                        B_new[i][j] = B_old[i][j];
                    }
                }
                for(i=0; i<row; i++){
                    for(j=0; j<K; j++){
                        A_new[i][j] = A_old[i][j];
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    grad = 0.0;
                    hess = 0.0;
                    for(k=0; k<col ; k++){
                        sum = exp(product[i][k]);
                        grad = grad + (sum - weight[i][k]) * B_new[j][k];
                        hess = hess + sum * B_new[j][k] * B_new[j][k];
                    }
                    if(hess != 0.0){
                        A_new[i][j] -= step_size * grad / hess;
                    }
                }
            }
            for(i=0; i<row; i++){
                sum = 0.0;
                for(j=0; j<K; j++){
                    sum = sum + A_new[i][j] * A_new[i][j];
                }
                if(sum > M*M){
                    for(j=0; j<K; j++){
                        A_new[i][j] = A_new[i][j] / sqrt(sum) * M;
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    grad = 0.0;
                    hess = 0.0;
                    for(k=0; k<row; k++){
                        sum = exp(product[k][j]);
                        grad += (sum - weight[k][j]) * A_new[k][i];
                        hess += sum * A_new[k][i] * A_new[k][i];
                    }
                    if(hess != 0.0){
                        B_new[i][j] -= step_size * grad / hess;
                    }
                }
            }
            for(j=0; j<col; j++){
                sum = 0.0;
                for(k=0; k<K; k++){
                    sum = sum + B_new[k][j] * B_new[k][j];
                }
                if(sum > M*M){
                    for(k=0; k<K; k++){
                        B_new[k][j] = B_new[k][j] / sqrt(sum) * M;
                    }
                }
            }
            f_new = 0.0;
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                    f_new = f_new - (weight[i][j] * product[i][j] - exp(product[i][j]));
                }
            }
            count_grad++;
            if(count_grad > max_grad){
                break;
            }
        }
        count_iter++;
        if(count_iter > max_iter){
            break;
        }
    }
    sum = 0.0;
    double product_true[row][col];
    int time_node;
    for(time_node=0; time_node<grid_num; time_node++){
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product_true[i][j] = 0.0;
                for(k=0; k<K_true; k++){
                    product_true[i][j] += rate_left_true(i,k) * rate_right_true(k,j);
                }
                if(time_node == 0 || time_node == grid_num - 1){
                    sum += 0.5 * (product_true[i][j] - product[i][j]) * (product_true[i][j] - product[i][j]);
                }else{
                    sum += (product_true[i][j] - product[i][j]) * (product_true[i][j] - product[i][j]);
                }
            }
        }
    }
    sum = sum * grid_size /row/col;
    List result = List(2);
    result[0] = f_new;
    result[1] = sum;
    return result;
} 
")
myFunction <- cppFunction("
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
//[[Rcpp::export]]
List poisson_linear(NumericMatrix event, int nrow, NumericMatrix T_mat, NumericMatrix rate_left_true, NumericMatrix rate_right_true, NumericMatrix coef, int row, int col, int K, int K_true, int max_grad, int max_iter, double step_size_initial, double M, double threshold, int grid_num, double t_min, double t_max) {
	  int i, j, k, count_iter, count_grad, row_event, col_event;
    double grad, hess, sum, f_new = 0.0, f_old = 0.0, step_size, grid_size;
    double product[row][col];
    double weight[row][col];
    double B_old[K][col], B_new[K][col], A_new[row][K], A_old[row][K];
    grid_size = (t_max - t_min)/(grid_num - 1);
    double A[row][K], B[K][col];
    if(K <= K_true){
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_new[i][j] = rate_left_true(i,j);
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = rate_right_true(i,j);
            }
        }
    }else{
        for(i=0; i<row; i++){
            for(j=0; j<K_true; j++){
                A[i][j] = rate_left_true(i,j);
            }
            for(j=K_true; j<K; j++){
                A[i][j] = 0.0;
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_new[i][j] = 0.0;
                for(k=0; k<K; k++){
                    A_new[i][j] += A[i][k] * T_mat(k,j);
                }
            }
        }
        for(i=0; i<K_true; i++){
            for(j=0; j<col; j++){
                B[i][j] = rate_right_true(i,j);
            }
        }
        for(i=K_true; i<K; i++){
            for(j=0; j<col; j++){
                B[i][j] = 0.0;
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = 0.0;
                for(k=0; k<K; k++){
                    B_new[i][j] += T_mat(k,i) * B[k][j];
                }
            }
        }
    }
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            weight[i][j] = 0.0;
        }
    }
    for(i=0;i<nrow;i++){
        row_event = event(i,0) - 1;
        col_event = event(i,1) - 1;
        weight[row_event][col_event] ++;
    }
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            product[i][j] = 0.0;
            for(k=0; k<K; k++){
                product[i][j] = product[i][j] + A_new[i][k] * B_new[k][j];
            }
            f_new = f_new - (weight[i][j] * product[i][j] - exp(product[i][j]));
        }
    }
    count_iter = 0;
    while(count_iter == 0|| fabs(f_new - f_old) > threshold){
        f_old = f_new;
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_old[i][j] = B_new[i][j];
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_old[i][j] = A_new[i][j];
            }
        }
        count_grad = 0;
        step_size = step_size_initial;
        while(f_new >= f_old){
            step_size = step_size * 0.5;
            if(count_grad > 0){
                for(i=0; i<K; i++){
                    for(j=0; j<col; j++){
                        B_new[i][j] = B_old[i][j];
                    }
                }
                for(i=0; i<row; i++){
                    for(j=0; j<K; j++){
                        A_new[i][j] = A_old[i][j];
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    grad = 0.0;
                    hess = 0.0;
                    for(k=0; k<col ; k++){
                        sum = exp(product[i][k]);
                        grad = grad + (sum - weight[i][k]) * B_new[j][k];
                        hess = hess + sum * B_new[j][k] * B_new[j][k];
                    }
                    if(hess != 0.0){
                        A_new[i][j] -= step_size * grad / hess;
                    }
                }
            }
            for(i=0; i<row; i++){
                sum = 0.0;
                for(j=0; j<K; j++){
                    sum = sum + A_new[i][j] * A_new[i][j];
                }
                if(sum > M*M){
                    for(j=0; j<K; j++){
                        A_new[i][j] = A_new[i][j] / sqrt(sum) * M;
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    grad = 0.0;
                    hess = 0.0;
                    for(k=0; k<row; k++){
                        sum = exp(product[k][j]);
                        grad += (sum - weight[k][j]) * A_new[k][i];
                        hess += sum * A_new[k][i] * A_new[k][i];
                    }
                    if(hess != 0.0){
                        B_new[i][j] -= step_size * grad / hess;
                    }
                }
            }
            for(j=0; j<col; j++){
                sum = 0.0;
                for(k=0; k<K; k++){
                    sum = sum + B_new[k][j] * B_new[k][j];
                }
                if(sum > M*M){
                    for(k=0; k<K; k++){
                        B_new[k][j] = B_new[k][j] / sqrt(sum) * M;
                    }
                }
            }
            f_new = 0.0;
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                    f_new = f_new - (weight[i][j] * product[i][j] - exp(product[i][j]));
                }
            }
            count_grad++;
            if(count_grad > max_grad){
                break;
            }
        }
        count_iter++;
        if(count_iter > max_iter){
            break;
        }
    }
    sum = 0.0;
    double product_true[row][col], t;
    int time_node;
    for(time_node=0; time_node<grid_num; time_node++){
        t = t_min + time_node * grid_size;
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product_true[i][j] = 0.0;
                for(k=0; k<K_true; k++){
                    product_true[i][j] += (rate_left_true(i,k) + t * coef(i,k))* rate_right_true(k,j);
                }
                if(time_node == 0 || time_node == grid_num - 1){
                    sum += 0.5 * (product_true[i][j] - product[i][j]) * (product_true[i][j] - product[i][j]);
                }else{
                    sum += (product_true[i][j] - product[i][j]) * (product_true[i][j] - product[i][j]);
                }
            }
        }
    }
    sum = sum * grid_size /row/col;
    List result = List(2);
    result[0] = f_new;
    result[1] = sum;
    return result;
} 
")
myFunction <- cppFunction("
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
//[[Rcpp::export]]
List poisson_periodic(NumericMatrix event, int nrow, NumericMatrix T_mat, NumericMatrix rate_left_true, NumericMatrix rate_right_true, NumericMatrix coef_scale, NumericMatrix coef_phase, double period, int row, int col, int K, int K_true, int max_grad, int max_iter, double step_size_initial, double M, double threshold, int grid_num, double t_min, double t_max) {
	  int i, j, k, count_iter, count_grad, row_event, col_event;
    double grad, hess, sum, f_new = 0.0, f_old = 0.0, step_size, grid_size;
    double product[row][col];
    double weight[row][col];
    double B_old[K][col], B_new[K][col], A_new[row][K], A_old[row][K];
    grid_size = (t_max - t_min)/(grid_num - 1);
    double A[row][K], B[K][col];
    if(K <= K_true){
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_new[i][j] = rate_left_true(i,j);
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = rate_right_true(i,j);
            }
        }
    }else{
        for(i=0; i<row; i++){
            for(j=0; j<K_true; j++){
                A[i][j] = rate_left_true(i,j);
            }
            for(j=K_true; j<K; j++){
                A[i][j] = 0.0;
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_new[i][j] = 0.0;
                for(k=0; k<K; k++){
                    A_new[i][j] += A[i][k] * T_mat(k,j);
                }
            }
        }
        for(i=0; i<K_true; i++){
            for(j=0; j<col; j++){
                B[i][j] = rate_right_true(i,j);
            }
        }
        for(i=K_true; i<K; i++){
            for(j=0; j<col; j++){
                B[i][j] = 0.0;
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = 0.0;
                for(k=0; k<K; k++){
                    B_new[i][j] += T_mat(k,i) * B[k][j];
                }
            }
        }
    }
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            weight[i][j] = 0.0;
        }
    }
    for(i=0;i<nrow;i++){
        row_event = event(i,0) - 1;
        col_event = event(i,1) - 1;
        weight[row_event][col_event] ++;
    }
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            product[i][j] = 0.0;
            for(k=0; k<K; k++){
                product[i][j] = product[i][j] + A_new[i][k] * B_new[k][j];
            }
            f_new = f_new - (weight[i][j] * product[i][j] - exp(product[i][j]));
        }
    }
    count_iter = 0;
    while(count_iter == 0|| fabs(f_new - f_old) > threshold){
        f_old = f_new;
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_old[i][j] = B_new[i][j];
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_old[i][j] = A_new[i][j];
            }
        }
        count_grad = 0;
        step_size = step_size_initial;
        while(f_new >= f_old){
            step_size = step_size * 0.5;
            if(count_grad > 0){
                for(i=0; i<K; i++){
                    for(j=0; j<col; j++){
                        B_new[i][j] = B_old[i][j];
                    }
                }
                for(i=0; i<row; i++){
                    for(j=0; j<K; j++){
                        A_new[i][j] = A_old[i][j];
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    grad = 0.0;
                    hess = 0.0;
                    for(k=0; k<col ; k++){
                        sum = exp(product[i][k]);
                        grad = grad + (sum - weight[i][k]) * B_new[j][k];
                        hess = hess + sum * B_new[j][k] * B_new[j][k];
                    }
                    if(hess != 0.0){
                        A_new[i][j] -= step_size * grad / hess;
                    }
                }
            }
            for(i=0; i<row; i++){
                sum = 0.0;
                for(j=0; j<K; j++){
                    sum = sum + A_new[i][j] * A_new[i][j];
                }
                if(sum > M*M){
                    for(j=0; j<K; j++){
                        A_new[i][j] = A_new[i][j] / sqrt(sum) * M;
                    }
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    grad = 0.0;
                    hess = 0.0;
                    for(k=0; k<row; k++){
                        sum = exp(product[k][j]);
                        grad += (sum - weight[k][j]) * A_new[k][i];
                        hess += sum * A_new[k][i] * A_new[k][i];
                    }
                    if(hess != 0.0){
                        B_new[i][j] -= step_size * grad / hess;
                    }
                }
            }
            for(j=0; j<col; j++){
                sum = 0.0;
                for(k=0; k<K; k++){
                    sum = sum + B_new[k][j] * B_new[k][j];
                }
                if(sum > M*M){
                    for(k=0; k<K; k++){
                        B_new[k][j] = B_new[k][j] / sqrt(sum) * M;
                    }
                }
            }
            f_new = 0.0;
            for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                    product[i][j] = 0.0;
                    for(k=0; k<K; k++){
                        product[i][j] += A_new[i][k] * B_new[k][j];
                    }
                    f_new = f_new - (weight[i][j] * product[i][j] - exp(product[i][j]));
                }
            }
            count_grad++;
            if(count_grad > max_grad){
                break;
            }
        }
        count_iter++;
        if(count_iter > max_iter){
            break;
        }
    }
    sum = 0.0;
    double product_true[row][col], t;
    int time_node;
    for(time_node=0; time_node<grid_num; time_node++){
        t = t_min + time_node * grid_size;
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product_true[i][j] = 0.0;
                for(k=0; k<K_true; k++){
                    product_true[i][j] += (rate_left_true(i,k) + coef_scale(i,k) * sin(2*M_PI/period*(t+coef_phase(i,k))))* rate_right_true(k,j);
                }
                if(time_node == 0 || time_node == grid_num - 1){
                    sum += 0.5 * (product_true[i][j] - product[i][j]) * (product_true[i][j] - product[i][j]);
                }else{
                    sum += (product_true[i][j] - product[i][j]) * (product_true[i][j] - product[i][j]);
                }
            }
        }
    }
    sum = sum * grid_size /row/col;
    List result = List(2);
    result[0] = f_new;
    result[1] = sum;
    return result;
} 
")

K_true <- 3
event_max <- 40000000
N_individual <- c(400,800,1200,1600,400,800,1200,1600)
N_item <- c(200,400,600,800,200,400,600,800)
n <- length(N_individual)/2
low <- 1.8
high <- 1.8
Mu_1 <- rep(-low, n*2)
Mu_2 <- rep(high, n*2)
Mu_3 <- rep(-high, n*2)
Mu_4 <- rep(low, n*2)
Mu_5 <- c(rep(-high, n), rep(-high/2, n))
Mu_6 <- c(rep(low, n), rep(low/2, n))
n_setting <- length(N_individual)
max_iter <- 200
max_grad <- 30
grid_num <- 30
threshold <- 10^(-6)
step_size_initial <- 1
M <- 6
period <- 1

Log_lik_int_poi <- array(rep(0,n_setting*3),dim = c(n_setting,3))
Loss_poi <- array(rep(0,n_setting*3),dim = c(n_setting,3))

#constant case
for(i in 1:n_setting){
  n_individual <- N_individual[i]
  n_item <- N_item[i]
  mu1 <- Mu_1[i]
  mu2 <- Mu_2[i]
  mu3 <- Mu_3[i]
  mu4 <- Mu_4[i]
  mu5 <- Mu_5[i]
  mu6 <- Mu_6[i]
  event <- matrix(rep(0,3*event_max),ncol=3)
  rate_log <- matrix(rep(0,n_individual*n_item),ncol=n_item)
  rate_log_left <- matrix(rep(0,K_true*n_individual),nrow=n_individual)
  rate_log_right <- matrix(rep(0,K_true*n_item),ncol=n_item)
  event_count <- generate_data_constant(event, n_individual, n_item, K_true, event_max, mu1, mu2, mu3, mu4, mu5, mu6, rate_log, rate_log_left, rate_log_right)
  transaction <- data.frame(household_key=event[(1:event_count),1],PRODUCT_ID=event[(1:event_count),2],time=event[(1:event_count),3])
  colnames(transaction) <- c("household_key","PRODUCT_ID","time")
  h <- 0.1*(min(c(n_item,n_individual))/(log(min(c(n_item,n_individual)))*log(min(c(n_item,n_individual)))))^(-1/5)
  t <- seq(h, 1-h, (1-2*h)/grid_num)
  K <- K_true
  T_mat <- randortho(K)
  result <- poisson_constant(as.matrix(transaction), nrow(transaction), T_mat, rate_log_left, rate_log_right, n_individual, n_item, K, K_true, max_grad, max_iter, step_size_initial, M, threshold, grid_num + 1, min(t), max(t))
  Log_lik_int_poi[i,1] <- result[[1]]
  Loss_poi[i,1] <- result[[2]]
  print(c(i,result[[2]]))
}

low <- 1.6
high <- 1.6
Mu_1 <- rep(-low, n*2)
Mu_2 <- rep(high, n*2)
Mu_3 <- rep(-high, n*2)
Mu_4 <- rep(low, n*2)
Mu_5 <- c(rep(-high, n), rep(-high/2, n))
Mu_6 <- c(rep(low, n), rep(low/2, n))

#linear case
for(i in 1:n_setting){
  n_individual <- N_individual[i]
  n_item <- N_item[i]
  mu1 <- Mu_1[i]
  mu2 <- Mu_2[i]
  mu3 <- Mu_3[i]
  mu4 <- Mu_4[i]
  mu5 <- Mu_5[i]
  mu6 <- Mu_6[i]
  event <- matrix(rep(0,3*event_max),ncol=3)
  rate_log <- matrix(rep(0,n_individual*n_item),ncol=n_item)
  rate_log_left <- matrix(rep(0,K_true*n_individual),nrow=n_individual)
  rate_log_right <- matrix(rep(0,K_true*n_item),ncol=n_item)
  coef <- matrix(runif((K_true)*n_individual,min=-4,max=4),nrow=n_individual)
  event_count <- generate_data_linear(event, n_individual, n_item, K_true, event_max, mu1, mu2, mu3, mu4, mu5, mu6, rate_log, rate_log_left, rate_log_right, coef)
  transaction <- data.frame(household_key=event[(1:event_count),1],PRODUCT_ID=event[(1:event_count),2],time=event[(1:event_count),3])
  colnames(transaction) <- c("household_key","PRODUCT_ID","time")
  h <- 0.1*(min(c(n_item,n_individual))/(log(min(c(n_item,n_individual)))*log(min(c(n_item,n_individual)))))^(-1/5)
  t <- seq(h, 1-h, (1-2*h)/grid_num)
  K <- K_true
  T_mat <- randortho(K)
  result <- poisson_linear(as.matrix(transaction), nrow(transaction), T_mat, rate_log_left, rate_log_right, coef, n_individual, n_item, K, K_true, max_grad, max_iter, step_size_initial, M, threshold, grid_num + 1, min(t), max(t))
  Log_lik_int_poi[i,2] <- result[[1]]
  Loss_poi[i,2] <- result[[2]]
  print(c(i,result[[2]]))
}

#periodic case
for(i in 1:n_setting){
  n_individual <- N_individual[i]
  n_item <- N_item[i]
  mu1 <- Mu_1[i]
  mu2 <- Mu_2[i]
  mu3 <- Mu_3[i]
  mu4 <- Mu_4[i]
  mu5 <- Mu_5[i]
  mu6 <- Mu_6[i]
  event <- matrix(rep(0,3*event_max),ncol=3)
  rate_log <- matrix(rep(0,n_individual*n_item),ncol=n_item)
  rate_log_left <- matrix(rep(0,K_true*n_individual),nrow=n_individual)
  rate_log_right <- matrix(rep(0,K_true*n_item),ncol=n_item)
  coef_scale <- matrix(runif(K_true*n_individual,min=-1.2,max=1.2),nrow=n_individual)
  coef_phase <- matrix(runif(K_true*n_individual,min=0,max=period),nrow=n_individual)
  event_count <- generate_data_periodic(event, n_individual, n_item, K_true, event_max, mu1, mu2, mu3, mu4, mu5, mu6, rate_log, rate_log_left, rate_log_right, coef_scale, coef_phase, period)
  transaction <- data.frame(household_key=event[(1:event_count),1],PRODUCT_ID=event[(1:event_count),2],time=event[(1:event_count),3])
  colnames(transaction) <- c("household_key","PRODUCT_ID","time")
  h <- 0.1*(min(c(n_item,n_individual))/(log(min(c(n_item,n_individual)))*log(min(c(n_item,n_individual)))))^(-1/5)
  t <- seq(h, 1-h, (1-2*h)/grid_num)
  K <- K_true
  T_mat <- randortho(K)
  result <- poisson_periodic(as.matrix(transaction), nrow(transaction), T_mat, rate_log_left, rate_log_right, coef_scale, coef_phase, period, n_individual, n_item, K, K_true, max_grad, max_iter, step_size_initial, M, threshold, grid_num + 1, min(t), max(t))
  Log_lik_int_poi[i,3] <- result[[1]]
  Loss_poi[i,3] <- result[[2]]
  print(c(i,result[[2]]))
}

rm(rate_log)
rm(event)
rm(transaction)
save.image(paste0("poisson_",id,".RData"))