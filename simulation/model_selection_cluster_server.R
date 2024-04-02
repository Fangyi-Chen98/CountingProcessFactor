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
int generate_data_constant_cluster(NumericMatrix event, int n_individual, int n_item, int K, int max_event, double mu1, double mu2, double mu3, double mu4, double mu5, double mu6, NumericMatrix Rate, NumericMatrix Rate_left, NumericMatrix Rate_right, int cluster){
    GetRNGstate();
    int i,j,m,k,event_count=0;
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
    int group, last;
    group = int(n_item / cluster);
    last = n_item - cluster * group;
    double rate_max, rate_current, rate_max_whole;
    for(i=0;i<n_individual;i++){
        for(j=0;j<group;j++){
            rate_max_whole = -10000.0;
            for(k=0; k<cluster; k++){
                rate_max = Rate(i,(j*cluster+k));
                if(rate_max > rate_max_whole){
                    rate_max_whole = rate_max;
                }
            }
            time_left = T;
            current_time = 0.0;
            while(time_left>0){
                double U = R::rexp(exp(-rate_max_whole));
                if(time_left > U){
                    current_time = current_time + U;
                    time_left = time_left - U;
                    for(k=0; k<cluster; k++){
                        rate_current = Rate(i,(j*cluster+k));
                        double U_1 = R::runif(0,1);
                        if(U_1 < exp(rate_current-rate_max_whole)){
                            event(event_count,0) = i+1;
                            event(event_count,1) = (j*cluster+k)+1;
                            event(event_count,2) = current_time;
                            event_count++;
                        }
                        if(event_count == max_event){
                            return -1;
                        }
                    }
                }else{
                    break;
                }
            }
        }
        if(last > 0){
            rate_max_whole = -10000.0;
            for(k=0; k<last; k++){
                rate_max = Rate(i,(j*cluster+k));
                if(rate_max > rate_max_whole){
                    rate_max_whole = rate_max;
                }
            }
            time_left = T;
            current_time = 0.0;
            while(time_left>0){
                double U = R::rexp(exp(-rate_max_whole));
                if(time_left > U){
                    current_time = current_time + U;
                    time_left = time_left - U;
                    for(k=0; k<last; k++){
                        rate_current = Rate(i,(j*cluster+k));
                        double U_1 = R::runif(0,1);
                        if(U_1 < exp(rate_current-rate_max_whole)){
                            event(event_count,0) = i+1;
                            event(event_count,1) = (j*cluster+k)+1;
                            event(event_count,2) = current_time;
                            event_count++;
                        }
                        if(event_count == max_event){
                            return -1;
                        }
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
int generate_data_linear_cluster(NumericMatrix event, int n_individual, int n_item, int K, int max_event, double mu1, double mu2, double mu3, double mu4, double mu5, double mu6, NumericMatrix Rate, NumericMatrix Rate_left, NumericMatrix Rate_right, NumericMatrix coef, int cluster){
    GetRNGstate();
    int i,j,m,k,event_count=0;
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
    int group, last, interval;
    group = int(n_item / cluster);
    last = n_item - cluster * group;
    double rate_max, rate_current, rate_max_whole, Rate_current[n_individual][n_item];
    int interval_num = 50;
    for(interval = 0; interval < interval_num; interval ++){
        for(i=0; i<n_individual; i++){
            for(j=0; j<n_item; j++){
                Rate_current[i][j] = Rate(i,j);
                for(k = 0; k < K; k++){
                    Rate_current[i][j] += interval*1.00/interval_num * coef(i,k)* Rate_right(k,j);
                }
            }
        }
        for(i=0;i<n_individual;i++){
            for(j=0;j<group;j++){
                rate_max_whole = -10000.0;
                for(k=0; k<cluster; k++){
                    rate_max = 0;
                    for(m=0; m<K; m++){
                        rate_max += coef(i,m)*X_right[m][(j*cluster+k)];
                    }
                    rate_max = fabs(rate_max)/interval_num + Rate_current[i][(j*cluster+k)];
                    if(rate_max > rate_max_whole){
                        rate_max_whole = rate_max;
                    }
                }
                time_left = T/interval_num;
                current_time = interval*1.00/interval_num;
                while(time_left>0){
                    double U = R::rexp(exp(-rate_max_whole));
                    if(time_left > U){
                        current_time = current_time + U;
                        time_left = time_left - U;
                        for(k=0; k<cluster; k++){
                            rate_current = Rate(i,(j*cluster+k));
                            for(m=0; m<K; m++){
                                rate_current += current_time * coef(i,m)*X_right[m][(j*cluster+k)];
                            }
                            double U_1 = R::runif(0,1);
                            if(U_1 < exp(rate_current-rate_max_whole)){
                                event(event_count,0) = i+1;
                                event(event_count,1) = (j*cluster+k)+1;
                                event(event_count,2) = current_time;
                                event_count++;
                            }
                            if(event_count == max_event){
                                return -1;
                            }
                        }
                    }else{
                        break;
                    }
                }
            }
            if(last > 0){
                rate_max_whole = -10000.0;
                for(k=0; k<last; k++){
                    rate_max = 0;
                    for(m=0; m<K; m++){
                        rate_max += coef(i,m)*X_right[m][(j*cluster+k)];
                    }
                    rate_max = fabs(rate_max)/interval_num + Rate_current[i][(j*cluster+k)];
                    if(rate_max > rate_max_whole){
                        rate_max_whole = rate_max;
                    }
                }
                time_left = T/interval_num;
                current_time = interval*1.00/interval_num;
                while(time_left>0){
                    double U = R::rexp(exp(-rate_max_whole));
                    if(time_left > U){
                        current_time = current_time + U;
                        time_left = time_left - U;
                        for(k=0; k<last; k++){
                            rate_current = Rate(i,(j*cluster+k));
                            for(m=0; m<K; m++){
                                rate_current += current_time * coef(i,m)*X_right[m][(j*cluster+k)];
                            }
                            double U_1 = R::runif(0,1);
                            if(U_1 < exp(rate_current-rate_max_whole)){
                                event(event_count,0) = i+1;
                                event(event_count,1) = (j*cluster+k)+1;
                                event(event_count,2) = current_time;
                                event_count++;
                            }
                            if(event_count == max_event){
                                return -1;
                            }
                        }
                    }else{
                        break;
                    }
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
int generate_data_periodic_cluster(NumericMatrix event, int n_individual, int n_item, int K, int max_event, double mu1, double mu2, double mu3, double mu4, double mu5, double mu6, NumericMatrix Rate, NumericMatrix Rate_left,NumericMatrix Rate_right, NumericMatrix coef_scale, NumericMatrix coef_phase, double period, int cluster){
    GetRNGstate();
    int i,j,m,k,event_count=0;
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
    int group, last, interval;
    group = int(n_item / cluster);
    last = n_item - cluster * group;
    double rate_max, rate_current, rate_max_whole, Rate_current[n_individual][n_item];
    int interval_num = 50;
    for(interval = 0; interval < interval_num; interval ++){   
        for(i=0; i<n_individual; i++){
            for(j=0; j<n_item; j++){
                Rate_current[i][j] = Rate(i,j);
                for(k = 0; k < K; k++){
                    Rate_current[i][j] += (coef_scale(i,k)*X_right[k][j])*sin(2*M_PI*(interval*1.00/interval_num+coef_phase(i,k))/period);
                }  
            }
        }
        for(i=0;i<n_individual;i++){
            for(j=0;j<group;j++){
                rate_max_whole = -10000.0;
                for(k=0; k<cluster; k++){
                    rate_max = Rate_current[i][(j*cluster+k)];
                    for(m=0; m<K; m++){
                        rate_max += fabs(coef_scale(i,m)*X_right[m][(j*cluster+k)])*1.00/interval_num;
                    }
                    if(rate_max > rate_max_whole){
                        rate_max_whole = rate_max;
                    }
                }
                time_left = T/interval_num;
                current_time = interval*1.00/interval_num;
                while(time_left>0){
                    double U = R::rexp(exp(-rate_max_whole));
                    if(time_left > U){
                        current_time = current_time + U;
                        time_left = time_left - U;
                        for(k=0; k<cluster; k++){
                            rate_current = Rate(i,(j*cluster+k));
                            for(m=0; m<K; m++){
                                rate_current += (coef_scale(i,m)*X_right[m][(j*cluster+k)])*sin(2*M_PI*(current_time+coef_phase(i,m))/period);
                            }
                            double U_1 = R::runif(0,1);
                            if(U_1 < exp(rate_current-rate_max_whole)){
                                event(event_count,0) = i+1;
                                event(event_count,1) = (j*cluster+k)+1;
                                event(event_count,2) = current_time;
                                event_count++;
                            }
                            if(event_count == max_event){
                                return -1;
                            }
                        }
                    }else{
                        break;
                    }
                }
            }
            if(last > 0){
                rate_max_whole = -10000.0;
                for(k=0; k<last; k++){
                    rate_max = Rate_current[i][(j*cluster+k)];
                    for(m=0; m<K; m++){
                        rate_max += fabs(coef_scale(i,m)*X_right[m][(j*cluster+k)])/interval_num;
                    }
                    if(rate_max > rate_max_whole){
                        rate_max_whole = rate_max;
                    }
                }
                time_left = T/interval_num;
                current_time = interval*1.00/interval_num;
                while(time_left>0){
                    double U = R::rexp(exp(-rate_max_whole));
                    if(time_left > U){
                        current_time = current_time + U;
                        time_left = time_left - U;
                        for(k=0; k<last; k++){
                            rate_current = Rate(i,(j*cluster+k));
                            for(m=0; m<K; m++){
                                rate_current += (coef_scale(i,m)*X_right[m][(j*cluster+k)])*sin(2*M_PI*(current_time+coef_phase(i,m))/period);
                            }
                            double U_1 = R::runif(0,1);
                            if(U_1 < exp(rate_current-rate_max_whole)){
                                event(event_count,0) = i+1;
                                event(event_count,1) = (j*cluster+k)+1;
                                event(event_count,2) = current_time;
                                event_count++;
                            }
                            if(event_count == max_event){
                                return -1;
                            }
                        }
                    }else{
                        break;
                    }
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
List simulation_constant(NumericMatrix event, int nrow, NumericMatrix T_mat, NumericMatrix rate_left_true, NumericMatrix rate_right_true, int row, int col, int K, int K_true, int max_grad, int max_iter, double step_size_initial, double M, double threshold, int grid_num, double t_min, double t_max, double h) {
	int i, j, k, count_iter, count_grad, event_num, time_node, row_event, col_event;
    double grad, hess, sum, f_new = 0.0, f_old = 0.0, step_size, error = 0.0, grid_size, t;
    double A_new[grid_num][row][K], A_old[grid_num][row][K];
    double product[grid_num][row][col], weight[grid_num][row][col], product_exp[grid_num][row][col];
    double B_old[K][col], B_new[K][col], B_grad[K][col], B_hess[K][col];
    grid_size = (t_max - t_min)/(grid_num - 1);
    double A[row][K], B[K][col];
    if(K <= K_true){
        for(time_node=0; time_node<grid_num; time_node++){
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    A_new[time_node][i][j] = rate_left_true(i,j);
                }
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = rate_right_true(i,j);
            }
        }
    }else{
        for(time_node=0; time_node<grid_num; time_node++){
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
                    A_new[time_node][i][j] = 0.0;
                    for(k=0; k<K; k++){
                        A_new[time_node][i][j] += A[i][k] * T_mat(k,j);
                    }
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
    for(time_node=0; time_node<grid_num; time_node++){
        t = t_min + time_node * grid_size;
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                weight[time_node][i][j] = 0.0;
            }
        }
        for(event_num=0; event_num<nrow; event_num++){
            row_event = event(event_num,0) - 1;
            col_event = event(event_num,1) - 1;
            if(fabs(event(event_num,2)-t) <= h){
              weight[time_node][row_event][col_event] += 3.0/4.0/h*(1 - (event(event_num,2)-t)*(event(event_num,2)-t)/h/h);
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product[time_node][i][j] = 0.0;
                for(k=0; k<K; k++){
                    product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                }
                product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                if(time_node == 0 || time_node == grid_num - 1){
                    f_new = f_new - 0.5 * (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                }else{
                    f_new = f_new - (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                }
            }
        }
    }
    f_new = f_new * grid_size;
    count_iter = 0;
    while(count_iter == 0|| error > threshold){
        f_old = f_new;
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_old[i][j] = B_new[i][j];
            }
        }
        for(time_node=0; time_node<grid_num; time_node++){
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    A_old[time_node][i][j] = A_new[time_node][i][j];
                }
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
                for(time_node=0; time_node<grid_num; time_node++){
                    for(i=0; i<row; i++){
                        for(j=0; j<K; j++){
                            A_new[time_node][i][j] = A_old[time_node][i][j];
                        }
                    }
                }
                for(time_node=0; time_node<grid_num; time_node++){
                    for(i=0; i<row; i++){
                        for(j=0; j<col; j++){
                            product[time_node][i][j] = 0.0;
                            for(k=0; k<K; k++){
                                product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                            }
                            product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                        }
                    }
                }
            }
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<row; i++){
                    for(j=0; j<K; j++){
                        grad = 0.0;
                        hess = 0.0;
                        for(k=0; k<col ; k++){
                            grad = grad + (product_exp[time_node][i][k] - weight[time_node][i][k]) * B_new[j][k];
                            hess = hess + product_exp[time_node][i][k] * B_new[j][k] * B_new[j][k];
                        }
                        if(hess != 0.0){
                            A_new[time_node][i][j] -= step_size * grad / hess;
                        }
                    }
                }
                for(i=0; i<row; i++){
                    sum = 0.0;
                    for(j=0; j<K; j++){
                        sum = sum + A_new[time_node][i][j] * A_new[time_node][i][j];
                    }
                    if(sum > M*M){
                        for(j=0; j<K; j++){
                            A_new[time_node][i][j] = A_new[time_node][i][j] / sqrt(sum) * M;
                        }
                    }
                }
                for(i=0; i<row; i++){
                    for(j=0; j<col; j++){
                        product[time_node][i][j] = 0.0;
                        for(k=0; k<K; k++){
                            product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                        }
                        product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    B_hess[i][j] = 0.0;
                    B_grad[i][j] = 0.0;
                }
            }
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<K; i++){
                    for(j=0; j<col; j++){
                        for(k=0; k<row; k++){
                            if(time_node == 0 || time_node == grid_num - 1){
                                B_grad[i][j] += 0.5 * (product_exp[time_node][k][j] - weight[time_node][k][j]) * A_new[time_node][k][i];
                                B_hess[i][j] += 0.5 * product_exp[time_node][k][j] * A_new[time_node][k][i] * A_new[time_node][k][i];
                            }else{
                                B_grad[i][j] += (product_exp[time_node][k][j] - weight[time_node][k][j]) * A_new[time_node][k][i];
                                B_hess[i][j] += product_exp[time_node][k][j] * A_new[time_node][k][i] * A_new[time_node][k][i];
                            }
                        }
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    if(B_hess[i][j] != 0.0){
                        B_new[i][j] -= step_size * B_grad[i][j] / B_hess[i][j];
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
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<row; i++){
                    for(j=0; j<col; j++){
                        product[time_node][i][j] = 0.0;
                        for(k=0; k<K; k++){
                            product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                        }
                        product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                        if(time_node == 0 || time_node == grid_num - 1){
                            f_new = f_new - 0.5 * (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                        }else{
                            f_new = f_new - (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                        }
                    }
                }
            }
            f_new = f_new * grid_size;
            count_grad++;
            if(count_grad > max_grad){
                break;
            }
        }
        count_iter++;
        if(count_iter > max_iter){
            break;
        }
        error = 0.0;
        for(time_node=0; time_node<grid_num; time_node++){
            for(i=0; i<row; i++){
                sum = 0.0;
                for(j=0; j<col; j++){
                    for(k=0; k<K; k++){
                        sum += A_new[time_node][i][k] * B_new[k][j] - A_old[time_node][i][k] * B_old[k][j];
                    }
                    if(time_node == 0 || time_node == grid_num - 1){
                        error = error + 0.5 * sum * sum;
                    }else{
                        error = error + sum * sum;
                    }
                }
            }
        }
        error = error * grid_size /row/col;
    }
    List result = List(2);
    result[0] = f_new;
    sum = 0.0;
    double product_true[row][col];
    for(time_node=0; time_node<grid_num; time_node++){
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product_true[i][j] = 0.0;
                for(k=0; k<K_true; k++){
                    product_true[i][j] += rate_left_true(i,k) * rate_right_true(k,j);
                }
                if(time_node == 0 || time_node == grid_num - 1){
                    sum += 0.5 * (product_true[i][j] - product[time_node][i][j]) * (product_true[i][j] - product[time_node][i][j]);
                }else{
                    sum += (product_true[i][j] - product[time_node][i][j]) * (product_true[i][j] - product[time_node][i][j]);
                }
            }
        }
    }
    sum = sum * grid_size /row/col;
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
List simulation_linear(NumericMatrix event, int nrow, NumericMatrix T_mat, NumericMatrix rate_left_true, NumericMatrix rate_right_true, NumericMatrix coef, int row, int col, int K, int K_true, int max_grad, int max_iter, double step_size_initial, double M, double threshold, int grid_num, double t_min, double t_max, double h) {
	int i, j, k, count_iter, count_grad, event_num, time_node, row_event, col_event;
    double grad, hess, sum, f_new = 0.0, f_old = 0.0, step_size, error = 0.0, grid_size, t;
    double A_new[grid_num][row][K], A_old[grid_num][row][K];
    double product[grid_num][row][col], weight[grid_num][row][col], product_exp[grid_num][row][col];
    double B_old[K][col], B_new[K][col], B_grad[K][col], B_hess[K][col];
    grid_size = (t_max - t_min)/(grid_num - 1);
    double A[row][K], B[K][col];
    if(K <= K_true){
        for(time_node=0; time_node<grid_num; time_node++){
            t = t_min + time_node * grid_size;
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    A_new[time_node][i][j] = rate_left_true(i,j) + t * coef(i,j);
                }
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = rate_right_true(i,j);
            }
        }
    }else{
        for(time_node=0; time_node<grid_num; time_node++){
            t = t_min + time_node * grid_size;
            for(i=0; i<row; i++){
                for(j=0; j<K_true; j++){
                    A[i][j] = rate_left_true(i,j) + t * coef(i,j);
                }
                for(j=K_true; j<K; j++){
                    A[i][j] = 0.0;
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    A_new[time_node][i][j] = 0.0;
                    for(k=0; k<K; k++){
                        A_new[time_node][i][j] += A[i][k] * T_mat(k,j);
                    }
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
    for(time_node=0; time_node<grid_num; time_node++){
        t = t_min + time_node * grid_size;
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                weight[time_node][i][j] = 0.0;
            }
        }
        for(event_num=0; event_num<nrow; event_num++){
            row_event = event(event_num,0) - 1;
            col_event = event(event_num,1) - 1;
            if(fabs(event(event_num,2)-t) <= h){
              weight[time_node][row_event][col_event] += 3.0/4.0/h*(1 - (event(event_num,2)-t)*(event(event_num,2)-t)/h/h);
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product[time_node][i][j] = 0.0;
                for(k=0; k<K; k++){
                    product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                }
                product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                if(time_node == 0 || time_node == grid_num - 1){
                    f_new = f_new - 0.5 * (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                }else{
                    f_new = f_new - (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                }
            }
        }
    }
    f_new = f_new * grid_size;
    count_iter = 0;
    while(count_iter == 0|| error > threshold){
        f_old = f_new;
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_old[i][j] = B_new[i][j];
            }
        }
        for(time_node=0; time_node<grid_num; time_node++){
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    A_old[time_node][i][j] = A_new[time_node][i][j];
                }
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
                for(time_node=0; time_node<grid_num; time_node++){
                    for(i=0; i<row; i++){
                        for(j=0; j<K; j++){
                            A_new[time_node][i][j] = A_old[time_node][i][j];
                        }
                    }
                }
                for(time_node=0; time_node<grid_num; time_node++){
                    for(i=0; i<row; i++){
                        for(j=0; j<col; j++){
                            product[time_node][i][j] = 0.0;
                            for(k=0; k<K; k++){
                                product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                            }
                            product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                        }
                    }
                }
            }
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<row; i++){
                    for(j=0; j<K; j++){
                        grad = 0.0;
                        hess = 0.0;
                        for(k=0; k<col ; k++){
                            grad = grad + (product_exp[time_node][i][k] - weight[time_node][i][k]) * B_new[j][k];
                            hess = hess + product_exp[time_node][i][k] * B_new[j][k] * B_new[j][k];
                        }
                        if(hess != 0.0){
                            A_new[time_node][i][j] -= step_size * grad / hess;
                        }
                    }
                }
                for(i=0; i<row; i++){
                    sum = 0.0;
                    for(j=0; j<K; j++){
                        sum = sum + A_new[time_node][i][j] * A_new[time_node][i][j];
                    }
                    if(sum > M*M){
                        for(j=0; j<K; j++){
                            A_new[time_node][i][j] = A_new[time_node][i][j] / sqrt(sum) * M;
                        }
                    }
                }
                for(i=0; i<row; i++){
                    for(j=0; j<col; j++){
                        product[time_node][i][j] = 0.0;
                        for(k=0; k<K; k++){
                            product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                        }
                        product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    B_hess[i][j] = 0.0;
                    B_grad[i][j] = 0.0;
                }
            }
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<K; i++){
                    for(j=0; j<col; j++){
                        for(k=0; k<row; k++){
                            if(time_node == 0 || time_node == grid_num - 1){
                                B_grad[i][j] += 0.5 * (product_exp[time_node][k][j] - weight[time_node][k][j]) * A_new[time_node][k][i];
                                B_hess[i][j] += 0.5 * product_exp[time_node][k][j] * A_new[time_node][k][i] * A_new[time_node][k][i];
                            }else{
                                B_grad[i][j] += (product_exp[time_node][k][j] - weight[time_node][k][j]) * A_new[time_node][k][i];
                                B_hess[i][j] += product_exp[time_node][k][j] * A_new[time_node][k][i] * A_new[time_node][k][i];
                            }
                        }
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    if(B_hess[i][j] != 0.0){
                        B_new[i][j] -= step_size * B_grad[i][j] / B_hess[i][j];
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
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<row; i++){
                    for(j=0; j<col; j++){
                        product[time_node][i][j] = 0.0;
                        for(k=0; k<K; k++){
                            product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                        }
                        product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                        if(time_node == 0 || time_node == grid_num - 1){
                            f_new = f_new - 0.5 * (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                        }else{
                            f_new = f_new - (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                        }
                    }
                }
            }
            f_new = f_new * grid_size;
            count_grad++;
            if(count_grad > max_grad){
                break;
            }
        }
        count_iter++;
        if(count_iter > max_iter){
            break;
        }
        error = 0.0;
        for(time_node=0; time_node<grid_num; time_node++){
            for(i=0; i<row; i++){
                sum = 0.0;
                for(j=0; j<col; j++){
                    for(k=0; k<K; k++){
                        sum += A_new[time_node][i][k] * B_new[k][j] - A_old[time_node][i][k] * B_old[k][j];
                    }
                    if(time_node == 0 || time_node == grid_num - 1){
                        error = error + 0.5 * sum * sum;
                    }else{
                        error = error + sum * sum;
                    }
                }
            }
        }
        error = error * grid_size /row/col;
    }
    List result = List(2);
    result[0] = f_new;
    sum = 0.0;
    double product_true[row][col];
    for(time_node=0; time_node<grid_num; time_node++){
        t = t_min + time_node * grid_size;
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product_true[i][j] = 0.0;
                for(k=0; k<K_true; k++){
                    product_true[i][j] += (rate_left_true(i,k) + t * coef(i,k))* rate_right_true(k,j);
                }
                if(time_node == 0 || time_node == grid_num - 1){
                    sum += 0.5 * (product_true[i][j] - product[time_node][i][j]) * (product_true[i][j] - product[time_node][i][j]);
                }else{
                    sum += (product_true[i][j] - product[time_node][i][j]) * (product_true[i][j] - product[time_node][i][j]);
                }
            }
        }
    }
    sum = sum * grid_size /row/col;
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
List simulation_periodic(NumericMatrix event, int nrow, NumericMatrix T_mat, NumericMatrix rate_left_true, NumericMatrix rate_right_true, NumericMatrix coef_scale, NumericMatrix coef_phase, double period, int row, int col, int K, int K_true, int max_grad, int max_iter, double step_size_initial, double M, double threshold, int grid_num, double t_min, double t_max, double h) {
	int i, j, k, count_iter, count_grad, event_num, time_node, row_event, col_event;
    double grad, hess, sum, f_new = 0.0, f_old = 0.0, step_size, error = 0.0, grid_size, t;
    double A_new[grid_num][row][K], A_old[grid_num][row][K];
    double product[grid_num][row][col], weight[grid_num][row][col], product_exp[grid_num][row][col];
    double B_old[K][col], B_new[K][col], B_grad[K][col], B_hess[K][col];
    grid_size = (t_max - t_min)/(grid_num - 1);
    double A[row][K], B[K][col];
    if(K <= K_true){
        for(time_node=0; time_node<grid_num; time_node++){
            t = t_min + time_node * grid_size;
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    A_new[time_node][i][j] = rate_left_true(i,j) + coef_scale(i,j) * sin(2*M_PI/period*(t+coef_phase(i,j)));;
                }
            }
        }
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_new[i][j] = rate_right_true(i,j);
            }
        }
    }else{
        for(time_node=0; time_node<grid_num; time_node++){
            t = t_min + time_node * grid_size;
            for(i=0; i<row; i++){
                for(j=0; j<K_true; j++){
                    A[i][j] = rate_left_true(i,j) + coef_scale(i,j) * sin(2*M_PI/period*(t+coef_phase(i,j)));;
                }
                for(j=K_true; j<K; j++){
                    A[i][j] = 0.0;
                }
            }
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    A_new[time_node][i][j] = 0.0;
                    for(k=0; k<K; k++){
                        A_new[time_node][i][j] += A[i][k] * T_mat(k,j);
                    }
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
    for(time_node=0; time_node<grid_num; time_node++){
        t = t_min + time_node * grid_size;
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                weight[time_node][i][j] = 0.0;
            }
        }
        for(event_num=0; event_num<nrow; event_num++){
            row_event = event(event_num,0) - 1;
            col_event = event(event_num,1) - 1;
            if(fabs(event(event_num,2)-t) <= h){
              weight[time_node][row_event][col_event] += 3.0/4.0/h*(1 - (event(event_num,2)-t)*(event(event_num,2)-t)/h/h);
            }
        }
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product[time_node][i][j] = 0.0;
                for(k=0; k<K; k++){
                    product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                }
                product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                if(time_node == 0 || time_node == grid_num - 1){
                    f_new = f_new - 0.5 * (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                }else{
                    f_new = f_new - (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                }
            }
        }
    }
    f_new = f_new * grid_size;
    count_iter = 0;
    while(count_iter == 0|| error > threshold){
        f_old = f_new;
        for(i=0; i<K; i++){
            for(j=0; j<col; j++){
                B_old[i][j] = B_new[i][j];
            }
        }
        for(time_node=0; time_node<grid_num; time_node++){
            for(i=0; i<row; i++){
                for(j=0; j<K; j++){
                    A_old[time_node][i][j] = A_new[time_node][i][j];
                }
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
                for(time_node=0; time_node<grid_num; time_node++){
                    for(i=0; i<row; i++){
                        for(j=0; j<K; j++){
                            A_new[time_node][i][j] = A_old[time_node][i][j];
                        }
                    }
                }
                for(time_node=0; time_node<grid_num; time_node++){
                    for(i=0; i<row; i++){
                        for(j=0; j<col; j++){
                            product[time_node][i][j] = 0.0;
                            for(k=0; k<K; k++){
                                product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                            }
                            product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                        }
                    }
                }
            }
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<row; i++){
                    for(j=0; j<K; j++){
                        grad = 0.0;
                        hess = 0.0;
                        for(k=0; k<col ; k++){
                            grad = grad + (product_exp[time_node][i][k] - weight[time_node][i][k]) * B_new[j][k];
                            hess = hess + product_exp[time_node][i][k] * B_new[j][k] * B_new[j][k];
                        }
                        if(hess != 0.0){
                            A_new[time_node][i][j] -= step_size * grad / hess;
                        }
                    }
                }
                for(i=0; i<row; i++){
                    sum = 0.0;
                    for(j=0; j<K; j++){
                        sum = sum + A_new[time_node][i][j] * A_new[time_node][i][j];
                    }
                    if(sum > M*M){
                        for(j=0; j<K; j++){
                            A_new[time_node][i][j] = A_new[time_node][i][j] / sqrt(sum) * M;
                        }
                    }
                }
                for(i=0; i<row; i++){
                    for(j=0; j<col; j++){
                        product[time_node][i][j] = 0.0;
                        for(k=0; k<K; k++){
                            product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                        }
                        product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    B_hess[i][j] = 0.0;
                    B_grad[i][j] = 0.0;
                }
            }
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<K; i++){
                    for(j=0; j<col; j++){
                        for(k=0; k<row; k++){
                            if(time_node == 0 || time_node == grid_num - 1){
                                B_grad[i][j] += 0.5 * (product_exp[time_node][k][j] - weight[time_node][k][j]) * A_new[time_node][k][i];
                                B_hess[i][j] += 0.5 * product_exp[time_node][k][j] * A_new[time_node][k][i] * A_new[time_node][k][i];
                            }else{
                                B_grad[i][j] += (product_exp[time_node][k][j] - weight[time_node][k][j]) * A_new[time_node][k][i];
                                B_hess[i][j] += product_exp[time_node][k][j] * A_new[time_node][k][i] * A_new[time_node][k][i];
                            }
                        }
                    }
                }
            }
            for(i=0; i<K; i++){
                for(j=0; j<col; j++){
                    if(B_hess[i][j] != 0.0){
                        B_new[i][j] -= step_size * B_grad[i][j] / B_hess[i][j];
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
            for(time_node=0; time_node<grid_num; time_node++){
                for(i=0; i<row; i++){
                    for(j=0; j<col; j++){
                        product[time_node][i][j] = 0.0;
                        for(k=0; k<K; k++){
                            product[time_node][i][j] += A_new[time_node][i][k] * B_new[k][j];
                        }
                        product_exp[time_node][i][j] = exp(product[time_node][i][j]);
                        if(time_node == 0 || time_node == grid_num - 1){
                            f_new = f_new - 0.5 * (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                        }else{
                            f_new = f_new - (weight[time_node][i][j] * product[time_node][i][j] - product_exp[time_node][i][j]);
                        }
                    }
                }
            }
            f_new = f_new * grid_size;
            count_grad++;
            if(count_grad > max_grad){
                break;
            }
        }
        count_iter++;
        if(count_iter > max_iter){
            break;
        }
        error = 0.0;
        for(time_node=0; time_node<grid_num; time_node++){
            for(i=0; i<row; i++){
                sum = 0.0;
                for(j=0; j<col; j++){
                    for(k=0; k<K; k++){
                        sum += A_new[time_node][i][k] * B_new[k][j] - A_old[time_node][i][k] * B_old[k][j];
                    }
                    if(time_node == 0 || time_node == grid_num - 1){
                        error = error + 0.5 * sum * sum;
                    }else{
                        error = error + sum * sum;
                    }
                }
            }
        }
        error = error * grid_size /row/col;
    }
    List result = List(2);
    result[0] = f_new;
    sum = 0.0;
    double product_true[row][col];
    for(time_node=0; time_node<grid_num; time_node++){
        t = t_min + time_node * grid_size;
        for(i=0; i<row; i++){
            for(j=0; j<col; j++){
                product_true[i][j] = 0.0;
                for(k=0; k<K_true; k++){
                    product_true[i][j] += (rate_left_true(i,k) + coef_scale(i,k) * sin(2*M_PI/period*(t+coef_phase(i,k))))* rate_right_true(k,j);
                }
                if(time_node == 0 || time_node == grid_num - 1){
                    sum += 0.5 * (product_true[i][j] - product[time_node][i][j]) * (product_true[i][j] - product[time_node][i][j]);
                }else{
                    sum += (product_true[i][j] - product[time_node][i][j]) * (product_true[i][j] - product[time_node][i][j]);
                }
            }
        }
    }
    sum = sum * grid_size /row/col;
    result[1] = sum;
    return result;
}
")
K_true <- 3
K_max <- 4
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
Cluster <- c(6,7,8,9,6,7,8,9)
n_setting <- length(N_individual)
max_iter <- 200
max_grad <- 30
grid_num <- 30
threshold <- 10^(-6)
step_size_initial <- 2
M <- 6
period <- 1

Log_lik_int_cluster <- array(rep(0,K_max*n_setting*3),dim = c(n_setting,K_max,3))
Loss_cluster <- array(rep(0,K_max*n_setting*3),dim = c(n_setting,K_max,3))

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
  cluster <- Cluster[i]
  event <- matrix(rep(0,3*event_max),ncol=3)
  rate_log <- matrix(rep(0,n_individual*n_item),ncol=n_item)
  rate_log_left <- matrix(rep(0,K_true*n_individual),nrow=n_individual)
  rate_log_right <- matrix(rep(0,K_true*n_item),ncol=n_item)
  event_count <- generate_data_constant_cluster(event, n_individual, n_item, K_true, event_max, mu1, mu2, mu3, mu4, mu5, mu6, rate_log, rate_log_left, rate_log_right, cluster)
  transaction <- data.frame(household_key=event[(1:event_count),1],PRODUCT_ID=event[(1:event_count),2],time=event[(1:event_count),3])
  colnames(transaction) <- c("household_key","PRODUCT_ID","time")
  h <- 0.1*(min(c(n_item,n_individual))/cluster)^(-1/5+0.01)
  t <- seq(h, 1-h, (1-2*h)/grid_num)
  for(k in 2:K_max){
    K <- k
    a <- gramSchmidt(diag(rep(1,K))+0.1)
    T_mat <- a$Q
    result <- simulation_constant(as.matrix(transaction), nrow(transaction), T_mat, rate_log_left, rate_log_right, n_individual, n_item, K, K_true, max_grad, max_iter, step_size_initial, M, threshold, grid_num + 1, min(t), max(t), h)
    Log_lik_int_cluster[i,k,1] <- result[[1]]
    Loss_cluster[i,k,1] <- result[[2]]
    print(c(i,k,result[[1]],result[[2]]))
  }
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
  cluster <- Cluster[i]
  event <- matrix(rep(0,3*event_max),ncol=3)
  rate_log <- matrix(rep(0,n_individual*n_item),ncol=n_item)
  rate_log_left <- matrix(rep(0,K_true*n_individual),nrow=n_individual)
  rate_log_right <- matrix(rep(0,K_true*n_item),ncol=n_item)
  coef <- matrix(runif((K_true)*n_individual,min=-4,max=4),nrow=n_individual)
  event_count <- generate_data_linear_cluster(event, n_individual, n_item, K_true, event_max, mu1, mu2, mu3, mu4, mu5, mu6, rate_log, rate_log_left, rate_log_right, coef, cluster)
  transaction <- data.frame(household_key=event[(1:event_count),1],PRODUCT_ID=event[(1:event_count),2],time=event[(1:event_count),3])
  colnames(transaction) <- c("household_key","PRODUCT_ID","time")
  h <- 0.1*(min(c(n_item,n_individual))/cluster)^(-1/5+0.01)
  t <- seq(h, 1-h, (1-2*h)/grid_num)
  for(k in 2:K_max){
    K <- k
    a <- gramSchmidt(diag(rep(1,K))+0.1)
    T_mat <- a$Q
    result <- simulation_linear(as.matrix(transaction), nrow(transaction), T_mat, rate_log_left, rate_log_right, coef, n_individual, n_item, K, K_true, max_grad, max_iter, step_size_initial, M, threshold, grid_num + 1, min(t), max(t), h)
    Log_lik_int_cluster[i,k,2] <- result[[1]]
    Loss_cluster[i,k,2] <- result[[2]]
    print(c(i,k,result[[1]],result[[2]]))
  }
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
  cluster <- Cluster[i]
  event <- matrix(rep(0,3*event_max),ncol=3)
  rate_log <- matrix(rep(0,n_individual*n_item),ncol=n_item)
  rate_log_left <- matrix(rep(0,K_true*n_individual),nrow=n_individual)
  rate_log_right <- matrix(rep(0,K_true*n_item),ncol=n_item)
  coef_scale <- matrix(runif(K_true*n_individual,min=-1.2,max=1.2),nrow=n_individual)
  coef_phase <- matrix(runif(K_true*n_individual,min=0,max=period),nrow=n_individual)
  event_count <- generate_data_periodic_cluster(event, n_individual, n_item, K_true, event_max, mu1, mu2, mu3, mu4, mu5, mu6, rate_log, rate_log_left, rate_log_right, coef_scale, coef_phase, period, cluster)
  transaction <- data.frame(household_key=event[(1:event_count),1],PRODUCT_ID=event[(1:event_count),2],time=event[(1:event_count),3])
  colnames(transaction) <- c("household_key","PRODUCT_ID","time")
  h <- 0.1*(min(c(n_item,n_individual))/cluster)^(-1/5+0.01)
  t <- seq(h, 1-h, (1-2*h)/grid_num)
  for(k in 2:K_max){
    K <- k
    a <- gramSchmidt(diag(rep(1,K))+0.1)
    T_mat <- a$Q
    result <- simulation_periodic(as.matrix(transaction), nrow(transaction), T_mat, rate_log_left, rate_log_right, coef_scale, coef_phase, period, n_individual, n_item, K, K_true, max_grad, max_iter, step_size_initial, M, threshold, grid_num + 1, min(t), max(t), h)
    Log_lik_int_cluster[i,k,3] <- result[[1]]
    Loss_cluster[i,k,3] <- result[[2]]
    print(c(i,k,result[[1]],result[[2]]))
  }
}

rm(rate_log)
rm(event)
rm(transaction)
save.image(paste0("model_selection_cluster_",id,".RData"))