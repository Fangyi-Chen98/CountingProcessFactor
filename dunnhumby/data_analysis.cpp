#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;
//[[Rcpp::export]]
NumericMatrix copy(NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericMatrix y(nrow, ncol);
    int i,j;
    for(i=0;i<nrow;i++){
        for(j=0;j<ncol;j++){
        y(i,j) = x(i,j);
        }
    }
    return y;
}
//[[Rcpp::export]]
List data_analysis(NumericMatrix event, int nrow, List initial_mat, int row, int col, int K, int max_grad, int max_iter, double step_size_initial, double M, double threshold, int grid_num, double t_min, double t_max, double h) {
	int i, j, k, count_iter, count_grad, event_num, time_node, row_event, col_event;
    double grad, hess, sum, f_new = 0.0, f_old = 0.0, step_size, grid_size, t;
    double A_new[grid_num][row][K], A_old[grid_num][row][K];
    double product[grid_num][row][col], weight[grid_num][row][col];
    double B_old[K][col], B_new[K][col], B_grad[K][col], B_hess[K][col];
    grid_size = (t_max - t_min)/(grid_num - 1);
    NumericMatrix A = copy(as<NumericMatrix>(initial_mat[0]));
    for(time_node=0; time_node<grid_num; time_node++){
        A = copy(as<NumericMatrix>(initial_mat[time_node]));
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A_new[time_node][i][j] = A(i,j);
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
    }
    NumericMatrix B = copy(as<NumericMatrix>(initial_mat[grid_num]));
    for(i=0; i<K; i++){
        for(j=0; j<col; j++){
            B_new[i][j] = B(i,j);
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
                if(time_node == 0 || time_node == grid_num - 1){
                    f_new = f_new - 0.5 * (weight[time_node][i][j] * product[time_node][i][j] - exp(product[time_node][i][j]));
                }else{
                    f_new = f_new - (weight[time_node][i][j] * product[time_node][i][j] - exp(product[time_node][i][j]));
                }
            }
        }
    }
    f_new = f_new * grid_size;
    count_iter = 0;
    printf("%f\n",f_new);
    while(count_iter == 0 || fabs(f_new - f_old) > threshold){
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
                            sum = exp(product[time_node][i][k]);
                            grad = grad + (sum - weight[time_node][i][k]) * B_new[j][k];
                            hess = hess + sum * B_new[j][k] * B_new[j][k];
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
                            sum = exp(product[time_node][k][j]);
                            if(time_node == 0 || time_node == grid_num - 1){
                                B_grad[i][j] += 0.5 * (sum - weight[time_node][k][j]) * A_old[time_node][k][i];
                                B_hess[i][j] += 0.5 * sum * A_old[time_node][k][i] * A_old[time_node][k][i];
                            }else{
                                B_grad[i][j] += (sum - weight[time_node][k][j]) * A_old[time_node][k][i];
                                B_hess[i][j] += sum * A_old[time_node][k][i] * A_old[time_node][k][i];
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
                        if(time_node == 0 || time_node == grid_num - 1){
                            f_new = f_new - 0.5 * (weight[time_node][i][j] * product[time_node][i][j] - exp(product[time_node][i][j]));
                        }else{
                            f_new = f_new - (weight[time_node][i][j] * product[time_node][i][j] - exp(product[time_node][i][j]));
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
        printf("%f %f\n",f_new, step_size);
    }
    List result = List(grid_num + 2);
    for(time_node=0; time_node<grid_num; time_node++){
        for(i=0; i<row; i++){
            for(j=0; j<K; j++){
                A(i,j) = A_new[time_node][i][j];
            }
        }
        result[time_node] = copy(as<NumericMatrix>(A));
    }
    for(i=0; i<K; i++){
        for(j=0; j<col; j++){
            B(i,j) = B_new[i][j];
        }
    }
    result[grid_num] = copy(as<NumericMatrix>(B));
    result[grid_num + 1] = f_new;
    return result;
}