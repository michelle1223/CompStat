# Computational Statistics HW4
# Chapter 5: Question 5.1, 5.3, 5.4
library("Rcpp")
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
sourceCpp("Computing_HW4_2020324011_민동은.cpp")

## Question 5.3 (a)
x = c(6.52, 8.32, 0.31, 2.82, 9.96, 0.14, 9.64)
xbar = mean(x)
k <- 7.84654
trapezoid(xbar, -100, 100, 10000)
k * trapezoid(xbar, -100, 100, 10000)

## Question 5.3 (b)
k * riemann(xbar, 2, 8, 500)
k * trapezoid(xbar, 2, 8, 500)
k * simpson(xbar, 2, 8, 500)

## Question 5.3 (c)
### First: ignore the singularity at 1
lower = exp(3)/(1+exp(3))
upper = 1
k * criemann(xbar, lower, upper, 500)
k * csimpson(xbar, lower, upper, 500)
### Second: fix the singularity at 1 by fixing the upper bound as a number slightly smaller than 1
upper = 1 - 0.00001
k * csimpson(xbar, lower, upper, 500)

## Question 5.3 (d)
lower = 10^(-10)
upper = 1/3
k * dsimpson(xbar, lower, upper, 1000)

## Question 5.4
a = 5
lower = (a-1)/a
upper = a-1
t_i0 = c(1:7)
for(i in 1:7){
  t_i0[i] = gtrapezoid(lower, upper, 2^(i-1))
}
t_i1 = c(1:6)
for(i in 1:6){
  t_i1[i] = (4*t_i0[i+1] - t_i0[i]) / (4-1)
}
t_i2 = c(1:5)
for(i in 1:5){
  t_i2[i] = (4*4*t_i1[i+1] - t_i1[i]) / (4^2-1)
}
t_i3 = c(1:4)
for(i in 1:4){
  t_i3[i] = ((4^3)*t_i2[i+1] - t_i2[i]) / (4^3-1)
}
t_i4 = c(1:3)
for(i in 1:3){
  t_i4[i] = ((4^4)*t_i3[i+1] - t_i3[i]) / (4^4-1)
}
t_i5 = c(1:2)
for(i in 1:2){
  t_i5[i] = ((4^5)*t_i4[i+1] - t_i4[i]) / (4^5-1)
}
t_66 = ((4^6)*t_i5[2] - t_i5[1]) / (4^6-1)
log(a)
