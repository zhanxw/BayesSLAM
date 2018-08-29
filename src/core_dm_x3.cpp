#include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include "common.h"

// [[Rcpp::export]]
Rcpp::List dm_model_estimator(IntegerMatrix Y, IntegerVector z, int iter, IntegerMatrix S, bool aggregate, bool fold_change, NumericMatrix b, NumericVector b_0) {
  // Read basic information from the data
  int n = Y.nrow();
  int p = Y.ncol();
  int K = max(z) + 1;
  int p_ext = S.nrow();
  int pp = S.ncol();
  if (!aggregate) 
  {
    p_ext = 0;
    pp = 0;
  }
  
  // Set hyperparameters
  double a_omega = 1.0;
  double b_omega = 1.0; 
  
  // For simulated data
  double a_0 = 2.0;
  double a = 2.0;
  double h_0 = 10;
  double h = 10;
  
  // Set algorithm settings
  int burn = iter*0.5;
  int E = (p_ext + p)*0.05;
  double tau = 1.0;
  double loga_lwr = -10.0;
  double loga_upp = 10.0;
  
  // Set temporary variables
  int it, i, ii, j, jj, k, e, gamma_temp, count_temp, count_temp_2;
  int count = 0;
  double hastings, temp, temp_2, alpha_temp;
  double accept_A = 0;
  double accept_gamma = 0;
  IntegerVector Y_sum(n);
  NumericVector A_sum(n);
  IntegerVector nk(K);
  IntegerVector gamma(p + p_ext);
  NumericVector gamma_ppi(p + p_ext);
  double gamma_sum = 0;
  IntegerVector gamma_sum_store(iter);
  NumericMatrix A(n, p + p_ext);
  NumericMatrix A_mean(n, p + p_ext);
  NumericMatrix logafc(iter - burn, p + p_ext);
  
  // Initialization
  for(k = 0; k < K; k++)
  {
    nk(k) = 0;
  }
  for(j = 0; j < p; j++)
  {
    gamma_ppi(j) = 0;
    gamma(j) = rbinom(1, 1, 0.05)(0);
    if(gamma(j) == 1)
    {
      gamma_sum++;
    }
    for(i = 0; i < n; i++)
    {
      A(i, j) = 1;
      A_mean(i, j) = 0;
      Y_sum(i) = Y_sum(i) + Y(i, j);
      A_sum(i) = A_sum(i) + A(i, j);
      if(j == 0)
      {
        nk(z(i)) = nk(z(i)) + 1;
      }
    }
  }
  if (aggregate)
  {
    for(j = 0; j < p_ext; j++)
    {
      gamma_ppi(p + j) = 0;
      gamma(p + j) = rbinom(1, 1, 0.05)(0);
      if(gamma(p + j) == 1)
      {
        gamma_sum++;
      }
      for(i = 0; i < n; i++)
      {
        A(i, p + j) = 0;
        A_mean(i, p + j) = 0;
        for(jj = 0; jj < pp; jj++)
        {
          if(S(j, jj) == 0)
          {
            break;
          }
          else
          {
            A(i, p + j) = A(i, p + j) + A(i, S(j, jj) - 1);
          }
        }
      }
    }
  }
  
  // MCMC
  for(it = 0; it < iter; it++)
  {
    // Update A
    for(i = 0; i < n; i++)
    {
      for(j = 0; j < p; j++)
      {
        //if(Y(i, j) != 0) 
        {
          //alpha_temp = exp(rnorm(1, log(A(i, j)), tau)(0));
          alpha_temp = exp(rnorm_trunc(log(A(i, j)), tau, loga_lwr, loga_upp));
          hastings = lgamma(A_sum(i) - A(i, j) + alpha_temp) - lgamma(Y_sum(i) + A_sum(i) - A(i, j) + alpha_temp) - lgamma(A_sum(i)) + lgamma(Y_sum(i) + A_sum(i));
          hastings = hastings + lgamma(Y(i, j) + alpha_temp) - lgamma(alpha_temp) - lgamma(Y(i, j) + A(i, j)) + lgamma(A(i, j));
          if(gamma(j) == 0)
          {
            temp = log(alpha_temp)*log(alpha_temp);
            temp_2 = log(alpha_temp);
            for(ii = 0; ii < n; ii++)
            {
              if(ii != i)
              {
                temp = temp + log(A(ii, j))*log(A(ii, j));
                temp_2 = temp_2 + log(A(ii, j));
              }
            }
            temp = temp - temp_2*temp_2/(n + 1/h_0);
            hastings = hastings - (a_0 + n*0.5)*log(b_0(j) + temp*0.5);
            temp = 0;
            temp_2 = 0;
            for(ii = 0; ii < n; ii++)
            {
              temp = temp + log(A(ii, j))*log(A(ii, j));
              temp_2 = temp_2 + log(A(ii, j));
            }
            temp = temp - temp_2*temp_2/(n + 1/h_0);
            hastings = hastings + (a_0 + n*0.5)*log(b_0(j) + temp*0.5);
          }
          else
          {
            temp = log(alpha_temp)*log(alpha_temp);
            temp_2 = log(alpha_temp);
            for(ii = 0; ii < n; ii++)
            {
              if(z[ii] == z[i] && ii != i)
              {
                temp = temp + log(A(ii, j))*log(A(ii, j));
                temp_2 = temp_2 + log(A(ii, j));
              }
            }
            temp = temp - temp_2*temp_2/(nk(z(i)) + 1/h);
            hastings = hastings - (a + nk(z(i))*0.5)*log(b(z(i), j) + temp*0.5);
            temp = 0;
            temp_2 = 0;
            for(ii = 0; ii < n; ii++)
            {
              if(z[ii] == z[i])
              {
                temp = temp + log(A(ii, j))*log(A(ii, j));
                temp_2 = temp_2 + log(A(ii, j));
              }
            }
            temp = temp - temp_2*temp_2/(nk(z(i)) + 1/h);
            hastings = hastings + (a + nk(z(i))*0.5)*log(b(z(i), j) + temp*0.5);
          }
          if (hastings >= log(double(rand()%10001)/10000))
          {
            A_sum(i) = A_sum(i) - A(i, j) + alpha_temp;
            A(i, j) = alpha_temp;
            if (it > burn) {
              accept_A++;
            }
          }
        }
      }
    }
    if (aggregate)
    {
      for(i = 0; i < n; i++)
      {
        for(j = 0; j < p_ext; j++)
        {
          A(i, p + j) = 0;
          for(jj = 0; jj < pp; jj++)
          {
            if(S(j, jj) == 0)
            {
              break;
            }
            else
            {
              A(i, p + j) = A(i, p + j) + A(i, S(j, jj) - 1);
            }
          }
        }
      }
    }
    
    
    
    // Update gamma
    for(e = 0; e < E; e++)
    {
      j = rand()%(p + p_ext);
      gamma_temp = 1 - gamma(j);
      if(gamma_temp == 0) // Delete
      {
        hastings = log(b_omega) - log(a_omega);
        temp = 0;
        temp_2 = 0;
        for(ii = 0; ii < n; ii++)
        {
          temp = temp + log(A(ii, j))*log(A(ii, j));
          temp_2 = temp_2 + log(A(ii, j));
        }
        temp = temp - temp_2*temp_2/(n + 1/h_0);
        hastings = hastings + (-log(n*h_0 + 1)*0.5 + lgamma(a_0 + n*0.5) - lgamma(a_0) + a_0*log(b_0(j)) - (a_0 + n*0.5)*log(b_0(j) + temp*0.5));
        for(k = 0; k < K; k++)
        {
          temp = 0;
          temp_2 = 0;
          for(ii = 0; ii < n; ii++)
          {
            if(z(ii) == k)
            {
              temp = temp + log(A(ii, j))*log(A(ii, j));
              temp_2 = temp_2 + log(A(ii, j));
            }
          }
          temp = temp - temp_2*temp_2/(nk(k) + 1/h);
          hastings = hastings - (-log(nk(k)*h + 1)*0.5 + lgamma(a + nk(k)*0.5) - lgamma(a) + a*log(b(k, j)) - (a + nk(k)*0.5)*log(b(k, j) + temp*0.5));
        }
      }
      else // Add
      {
        hastings = log(a_omega) - log(b_omega);
        for(k = 0; k < K; k++)
        {
          temp = 0;
          temp_2 = 0;
          for(ii = 0; ii < n; ii++)
          {
            if(z(ii) == k)
            {
              temp = temp + log(A(ii, j))*log(A(ii, j));
              temp_2 = temp_2 + log(A(ii, j));
            }
          }
          temp = temp - temp_2*temp_2/(nk(k) + 1/h);
          hastings = hastings + (-log(nk(k)*h + 1)*0.5 + lgamma(a + nk(k)*0.5) - lgamma(a) + a*log(b(k, j)) - (a + nk(k)*0.5)*log(b(k, j) + temp*0.5));
        }
        temp = 0;
        temp_2 = 0;
        for(ii = 0; ii < n; ii++)
        {
          temp = temp + log(A(ii, j))*log(A(ii, j));
          temp_2 = temp_2 + log(A(ii, j));
        }
        temp = temp - temp_2*temp_2/(n + 1/h_0);
        hastings = hastings - (-log(n*h_0 + 1)*0.5 + lgamma(a_0 + n*0.5) - lgamma(a_0) + a_0*log(b_0(j)) - (a_0 + n*0.5)*log(b_0(j) + temp*0.5));
      }
      
      if (hastings >= log(double(rand()%10001)/10000))
      {
        gamma(j) = gamma_temp;
        if(gamma_temp == 1) // Add
        {
          gamma_sum++;
        }
        else // Delete
        {
          gamma_sum--;
        }
        if(it > burn) {
          accept_gamma++;
        }
      }
    }
    
    // Monitor the process
    if(it*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    gamma_sum_store(it) = gamma_sum;
    if(it > burn) {
      for(j = 0; j < p + p_ext; j++)
      {
        gamma_ppi(j) = gamma_ppi(j) + gamma(j);
        temp = 0;
        temp_2 = 0;
        count_temp = 0;
        count_temp_2 = 0;
        for(i = 0; i < n; i++)
        {
          A_mean(i, j) = A_mean(i, j) + A(i, j);
          if(fold_change)
          {
            if(z(i) == 0)
            {
              temp = temp + A(i, j);
              count_temp++;
            }
            else if(z(i) == 1)
            {
              temp_2 = temp_2 + A(i, j);
              count_temp_2++;
            }
          }
        }
        if(fold_change)
        {
          logafc(it - burn - 1, j) = log(temp_2/count_temp_2) - log(temp/count_temp);
        }
      }
    }
  }
  
  accept_A = accept_A/n/p/(iter - burn);
  accept_gamma = accept_gamma/E/(iter - burn);
  for(j = 0; j < p + p_ext; j++)
  {
    gamma_ppi(j) = gamma_ppi(j)/(iter - burn);
    for(i = 0; i < n; i++)
    {
      if(j < p)
      {
        A_mean(i, j) = A_mean(i, j)/(iter - burn);
      }
    }
  }
  /*
  if (aggregate)
  {
    for(i = 0; i < n; i++)
    {
      for(j = 0; j < p_ext; j++)
      {
        A_mean(i, p + j) = 0;
        for(jj = 0; jj < pp; jj++)
        {
          if(S(j, jj) == 0)
          {
            break;
          }
          else
          {
            A_mean(i, p + j) = A_mean(i, p + j) + A_mean(i, S(j, jj) - 1);
          }
        }
      }
    }
  }
  */
  if(fold_change)
  {
    return Rcpp::List::create(Rcpp::Named("logafc") = logafc, Rcpp::Named("A_mean") = A_mean, Rcpp::Named("gamma_ppi") = gamma_ppi, Rcpp::Named("gamma_sum") = gamma_sum_store, Rcpp::Named("accept_A") = accept_A, Rcpp::Named("accept_gamma") = accept_gamma);
  }
  else
  {
    return Rcpp::List::create(Rcpp::Named("A_mean") = A_mean, Rcpp::Named("gamma_ppi") = gamma_ppi, Rcpp::Named("gamma_sum") = gamma_sum_store, Rcpp::Named("accept_A") = accept_A, Rcpp::Named("accept_gamma") = accept_gamma);
  }
}
