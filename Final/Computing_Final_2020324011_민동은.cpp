#include <Rcpp.h>
using namespace Rcpp;


// Question 3
// [[Rcpp::export]]
double lambda1(double x, double w, double a){
  LogicalVector e1 = (0.0 < x) && (x <= 1.5);
  LogicalVector e2 = (1.5 < x) && (x <= 3.0);
  LogicalVector e3 = (3.0 < x) && (x <= 4.5);
  LogicalVector e4 = (4.5 < x) && (x <= 6.0);
  NumericVector eta1 = 0.07116282*ifelse(e1, 1.0, 0.0) + 0.07766563*ifelse(e2, 1.0, 0.0) + 
    0.1052774*ifelse(e3, 1.0, 0.0) + 0.1061366*ifelse(e4, 1.0, 0.0);
  NumericVector l1 = eta1*exp(-0.1987954*w + -0.6744738*a + 0.1701579*w*a);
  return l1[0];
}

// [[Rcpp::export]]
double lambda2(double x, double w, double a){
  LogicalVector e1 = (0.0 < x) && (x <= 1.5);
  LogicalVector e2 = (1.5 < x) && (x <= 3.0);
  LogicalVector e3 = (3.0 < x) && (x <= 4.5);
  LogicalVector e4 = (4.5 < x) && (x <= 6.0);
  NumericVector eta2 = 0.02380466*ifelse(e1, 1.0, 0.0) + 0.02865413*ifelse(e2, 1.0, 0.0) + 
    0.03215047*ifelse(e3, 1.0, 0.0) + 0.03584044*ifelse(e4, 1.0, 0.0);
  NumericVector l2 = eta2*exp(-0.2518877*w + 0.3991342*a + 0.25158*w*a);
  return l2[0];
}

// [[Rcpp::export]]
double l1trapezoid(double W, double A, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sum = 0.0;
  for (i=1; i<n; i++){
    sum += lambda1(a+i*h, W, A);
  }
  return h*lambda1(a,W,A)/2 + h*sum + h*lambda1(b,W,A)/2;
}

// [[Rcpp::export]]
// n is an even number!!
double l1simpson(double W, double A, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sum = 0.0;
  for (i=1; i<(n/2)+1; i++){
    sum += lambda1(a+(2*i-2)*h, W, A) + 4*lambda1(a+(2*i-1)*h, W, A) + lambda1(a+(2*i*h), W, A);
  }
  return (h/3)*sum;
}

// [[Rcpp::export]]
double l2trapezoid(double W, double A, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sum = 0.0;
  for (i=1; i<n; i++){
    sum += lambda2(a+i*h, W, A);
  }
  return h*lambda2(a,W,A)/2 + h*sum + h*lambda2(b,W,A)/2;
}

// [[Rcpp::export]]
// n is an even number!!
double l2simpson(double W, double A, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sum = 0.0;
  for (i=1; i<(n/2)+1; i++){
    sum += lambda2(a+(2*i-2)*h, W, A) + 4*lambda2(a+(2*i-1)*h, W, A) + lambda2(a+(2*i*h), W, A);
  }
  return (h/3)*sum;
}

// [[Rcpp::export]]
double u1t(double X, double W, double A, int n){
  double u = exp(-l1trapezoid(W,A,0,X,n))*exp(-l2trapezoid(W,A,0,X,n))*lambda1(X,W,A);
  return u;
}

// [[Rcpp::export]]
double u2t(double X, double W, double A, int n){
  double u = exp(-l1trapezoid(W,A,0,X,n))*exp(-l2trapezoid(W,A,0,X,n))*lambda2(X,W,A);
  return u;
}

// [[Rcpp::export]]
double u1s(double X, double W, double A, int n){
  double u = exp(-l1simpson(W,A,0,X,n))*exp(-l2simpson(W,A,0,X,n))*lambda1(X,W,A);
  return u;
}

// [[Rcpp::export]]
double u2s(double X, double W, double A, int n){
  double u = exp(-l1simpson(W,A,0,X,n))*exp(-l2simpson(W,A,0,X,n))*lambda2(X,W,A);
  return u;
}

// [[Rcpp::export]]
double u1trapezoid(double X, double W, double A, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sum = 0.0;
  for (i=1; i<n; i++){
    sum += u1t(a+i*h, W, A, n);
  }
  return h*u1t(a,W,A,n)/2 + h*sum + h*u1t(b,W,A,n)/2;
}

// [[Rcpp::export]]
double u2trapezoid(double X, double W, double A, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sum = 0.0;
  for (i=1; i<n; i++){
    sum += u2t(a+i*h, W, A, n);
  }
  return h*u2t(a,W,A,n)/2 + h*sum + h*u2t(b,W,A,n)/2;
}

// [[Rcpp::export]]
double u1simpson(double X, double W, double A, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sum = 0.0;
  for (i=1; i<(n/2)+1; i++){
    sum += u1s(a+(2*i-2)*h,W,A,n) + 4*u1s(a+(2*i-1)*h,W,A,n) + u1s(a+(2*i*h),W,A,n);
  }
  return (h/3)*sum;
}

// [[Rcpp::export]]
double u2simpson(double X, double W, double A, double a, double b, int n){
  int i;
  double h = (b-a)/n;
  double sum = 0.0;
  for (i=1; i<(n/2)+1; i++){
    sum += u2s(a+(2*i-2)*h,W,A,n) + 4*u2s(a+(2*i-1)*h,W,A,n) + u2s(a+(2*i*h),W,A,n);
  }
  return (h/3)*sum;
}


// Question 4
// [[Rcpp::export]]
IntegerMatrix changematrix(IntegerMatrix table){
  int i; int j;
  int nrows = table.nrow(); 
  int ncols = table.ncol();
  IntegerMatrix mat(nrows+2, ncols+2);
  for (i=0; i<nrows; i++){
    for (j=0; j<ncols; j++){
      mat(i+1,j+1) = table(i,j);
    }
  }
  return mat;
}

// [[Rcpp::export]]
double findlength(IntegerMatrix table, int a){
  int n = 0;
  int nrows = table.nrow();
  int ncols = table.ncol();
  for (int i=0; i<nrows; i++){
    for (int j=0; j<ncols; j++){
      if (table(i,j)==a){
        n+=1;} } }
  return n;
}

// [[Rcpp::export]]
int neighbor(IntegerMatrix table, int i, int j){
  IntegerMatrix table1 = changematrix(table);
  int temp = table1(i+1,j) + table1(i+1,j+2) + table1(i,j+1) + table1(i+2,j+1);
  return temp;
}

// [[Rcpp::export]]
double boot(IntegerMatrix table, int i, int j, double a, double b){
  int x = table(i,j);
  if (x != 0){
    double s1 = findlength(table, 1)-findlength(table, -1);
    double s2 = s1*neighbor(table, i, j);
    double density1 = exp(a*s1 + (b/2)*s2);
    double s3 = findlength(table, 1)-findlength(table, -1);
    double s4 = s3*(-neighbor(table, i, j));
    double density2 = exp(a*s3 + (b/2)*s4);
    double density = density1 / (density1 + density2);
    double u = Rcpp::runif(1, 0.0, 1.0)[0];
    if (u < density){
      x = 1;
    } else{
      x = -1;
      }
  }
  return x;
}

// [[Rcpp::export]]
double test(IntegerMatrix table, int i, int j, double a, double b){
  double s1 = findlength(table, 1)-findlength(table, -1);
  double s2 = s1*neighbor(table, i, j);
  double density1 = exp(a*s1 + (b/2)*s2);
  double s3 = findlength(table, 1)-findlength(table, -1);
  double s4 = s3*(-neighbor(table, i, j));
  double density2 = exp(a*s3 + (b/2)*s4);
  double density = density1 / (density1 + density2);
  return density;
}

// [[Rcpp::export]]
IntegerMatrix bootstrap(IntegerMatrix table, double a, double b){
  int i; int j;
  int nrows = table.nrow();
  int ncols = table.ncol();
  for (i=0; i<nrows; i++){
    for (j=0; j<ncols; j++){
      table(i,j) = boot(table,i,j,a,b);
    }
  } return table;
}

// [[Rcpp::export]]
double rmse_s1(IntegerMatrix original, IntegerMatrix table){
  int s1 = findlength(original, 1) - findlength(original, -1);
  int s1_b = findlength(table, 1) - findlength(table, -1);
  double bias = std::pow((s1 - s1_b),2) / 2293;
  double var = std::pow((s1_b - (s1_b/2293)), 2) / 2293;
  return std::sqrt(bias+var);
}

// [[Rcpp::export]]
NumericVector rmse_s2(IntegerMatrix original, IntegerMatrix table){
  NumericVector x(2);
  int nrows = original.nrow();
  int ncols = original.ncol();
  double s2o = 0.0; double s2b = 0.0;
  for (int i=0; i<nrows; i++){
    for (int j=0; j<ncols; j++){
      double s2 = original(i,j)*neighbor(original,i,j);
      double s2_b = table(i,j)*neighbor(table,i,j);
      s2o += s2;
      s2b += s2_b; } }
  x[0] = s2o;
  x[1] = s2b;
  return x;
}

