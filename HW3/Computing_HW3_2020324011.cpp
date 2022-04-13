#include <Rcpp.h>
using namespace Rcpp;

// Example: Multinomial Distribution
// [[Rcpp::export]]
NumericVector multinom_e(NumericVector x, NumericVector p){
  NumericVector y(4);
  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2]*p[2]/(p[2]+p[3]);
  y[3] = x[2]*p[3]/(p[2]+p[3]);
  return y;
}

// [[Rcpp::export]]
double multinom_m(NumericVector x){
  double theta = (x[1]+x[2])/(x[0]+x[1]+x[2]);
  return theta;
}


// Example 4.2 Peppered Moths
// [[Rcpp::export]]
NumericVector allele_e(NumericVector x, NumericVector p){
  NumericVector n(6);
  n[0] = (x[0]*std::pow(p[0],2)) / (std::pow(p[0],2) + 2*p[0]*p[1] + 2*p[0]*p[2]);  // n_cc
  n[1] = (2*x[0]*p[0]*p[1]) / (std::pow(p[0],2) + 2*p[0]*p[1] + 2*p[0]*p[2]);  // n_ci
  n[2] = (2*x[0]*p[0]*p[2]) / (std::pow(p[0],2) + 2*p[0]*p[1] + 2*p[0]*p[2]);  // n_ct
  n[3] = (x[1]*std::pow(p[1],2)) / (std::pow(p[1],2) + 2*p[1]*p[2]);  // n_ii
  n[4] = (2*x[1]*p[1]*p[2]) / (std::pow(p[1],2) + 2*p[1]*p[2]);  // n_it
  n[5] = x[2];  // n_tt
  return n;
}

// [[Rcpp::export]]
NumericVector allele_m(NumericVector x, NumericVector n){
  NumericVector p(3);
  double sumx = 0.0;
  int i;
  int k = x.size();
  for (i=0; i<k; i++){
    sumx += x(i);
  }
  p[0] = (2*n[0]+n[1]+n[2])/(2*sumx);
  p[1] = (2*n[3]+n[4]+n[1])/(2*sumx);
  p[2] = (2*n[5]+n[2]+n[4])/(2*sumx);
  return p;
}


// Exercise 4.1 (a),(b)
// [[Rcpp::export]]
NumericVector moth_e(NumericVector x, NumericVector p){
  NumericVector n(6);
  n[0] = (x[0]*std::pow(p[0],2))/(std::pow(p[0],2) + 2*p[0]*p[1] + 2*p[0]*p[2]);  // n_cc
  n[1] = (2*x[0]*p[0]*p[1])/(std::pow(p[0],2) + 2*p[0]*p[1] + 2*p[0]*p[2]);  // n_ci
  n[2] = (2*x[0]*p[0]*p[2])/(std::pow(p[0],2) + 2*p[0]*p[1] + 2*p[0]*p[2]);  // n_ct
  n[3] = (x[1]*std::pow(p[1],2))/(std::pow(p[1],2) + 2*p[1]*p[2]) 
    + (x[3]*std::pow(p[1],2))/(std::pow(p[1],2) + 2*p[1]*p[2] + std::pow(p[2],2));  // n_ii
  n[4] = (2*x[1]*p[1]*p[2])/(std::pow(p[1],2) + 2*p[1]*p[2])
    + (x[3]*2*p[1]*p[2])/(std::pow(p[1],2) + 2*p[1]*p[2] + std::pow(p[2],2));  // n_it
  n[5] = x[2] + (x[3]*std::pow(p[2],2))/(std::pow(p[1],2) + 2*p[1]*p[2] + std::pow(p[2],2));  // n_tt
  return n;
}

// [[Rcpp::export]]
NumericVector moth_m(NumericVector x, NumericVector n){
  NumericVector p(3);
  double sumx = 0.0;
  int i;
  int k = x.size();
  for (i=0; i<k; i++){
    sumx += x(i);
  }
  p[0] = (2*n[0]+n[1]+n[2])/(2*sumx);
  p[1] = (2*n[3]+n[4]+n[1])/(2*sumx);
  p[2] = 1 - p[0] - p[1];
  return p;
}


// Exercise 4.2 (a),(b)
// x = (n_0, n_1, ... , n_16) and sum(x) = 1500, p = (alpha, beta, mu, lambda)
NumericVector concat(NumericVector a, NumericVector b){
  int n = b.size();
  for (int i=0; i<n; i++){
    a.push_back(b[i]);
  }
  return a;
}

// [[Rcpp::export]]
NumericVector hiv_e(NumericVector x, NumericVector p){
  int i;
  NumericVector pi_i(17);
  double n_z;
  NumericVector n_t(17);
  NumericVector n_p(17);
  pi_i[0] = p[0] + p[1]*exp(-p[2]) + (1-p[0]-p[1])*exp(-p[3]);
  for (i=1; i<17; i++){
    pi_i[i] = p[1]*exp(-p[2])*std::pow(p[2],i) + (1-p[0]-p[1])*exp(-p[3])*std::pow(p[3],i);
  }
  n_z = x[0]*p[0]/pi_i[0];
  for (i=0; i<17; i++){
    n_t[i] = x[i]*p[1]*std::pow(p[2],i)*exp(-p[2])/pi_i[i];
  }
  for (i=0; i<17; i++){
    n_p[i] = x[i]*(1-p[0]-p[1])*std::pow(p[3],i)*exp(-p[3])/pi_i[i];
  }
  NumericVector N = concat(n_t, n_p);
  N.push_front(n_z);
  return N;
}

// [[Rcpp::export]]
NumericVector hiv_m(NumericVector x){
  NumericVector p(4);  // p = (alpha, beta, mu, lambda)
  p[0] = x[0]/1500;
  int i;
  double sumt = 0.0, sumit = 0.0, sump = 0.0, sumip = 0.0;
  for (i=1; i<18; i++){
    sumt += x[i];
  }
  for (i=1; i<18; i++){
    sumit += (i-1)*x[i];
  }
  for (i=18; i<35; i++){
    sump += x[i];
  }
  for (i=18; i<35; i++){
    sumip += (i-18)*x[i];
  }
  p[1] = sumt/1500;
  p[2] = sumit/sumt;
  p[3] = sumip/sump;
  return p;
}

