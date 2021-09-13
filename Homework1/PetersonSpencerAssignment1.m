%% Problem 1
clear; clc; close all
x = linspace(0,2, 1001);
y = exp(x);

figure(101)
plot(x,y)
title('Simple Exponential Function')

%% Problem 2
clc
%close 201
% Evenly spaced points

all_test_cases = [6 11 21 41 81 161 321 641];

% figure(201)
% hold off
% plot(x,y);
% hold on;
k = 1;

norm_2_even = zeros(1,length(all_test_cases));
norm_inf_even = zeros(1,length(all_test_cases));

norm_2_cheby = zeros(1,length(all_test_cases));
norm_inf_cheby = zeros(1,length(all_test_cases));

vander_time = zeros(1,length(all_test_cases));
vander_time_cheby = zeros(1,length(all_test_cases));

for points = all_test_cases
   x_sampled = linspace(0,2,points);
   y_sampled = exp(x_sampled);
   
   tic
   V = vander(x_sampled);
   
   coeffs = V/y_sampled;
    
   
   result = polyval(coeffs,x);
   vander_time(k) = toc;
   
   norm_2_even(k) = norm(result-y,2);
   norm_inf_even(k) = norm(result-y, inf);
   
   n = points;
   x_chebyshev = 1 - cos(pi*(0:n-1)/(n-1));
   y_chebyshev = exp(x_chebyshev);
   
   tic
   V_cheby = vander(x_chebyshev);
   cheby_c = V_cheby/y_chebyshev;
   
   result_cheby = polyval(cheby_c,x);
   vander_time_cheby(k) = toc;
   
   norm_2_cheby(k) = norm(result_cheby - y,2);
   norm_inf_cheby(k) = norm(result_cheby - y,inf);
   
   k = k+1;
end

%% Problem 3


norm_2_lag = zeros(1,length(all_test_cases));
norm_inf_lag = zeros(1,length(all_test_cases));

norm_2_lag_cheby = zeros(1,length(all_test_cases));
norm_inf_lag_cheby = zeros(1,length(all_test_cases));


lag_time = zeros(1,length(all_test_cases));
lag_time_cheby = zeros(1,length(all_test_cases));

k = 1;
for points = all_test_cases
   x_lag = linspace(0,2,points);
   y_lag = exp(x_lag);
   
   tic
   result_lagrange = polyinterp(x_lag,y_lag, x);
   lag_time(k) = toc;
   
   norm_2_lag(k) = norm(result_lagrange-y,2);
   norm_inf_lag(k) = norm(result_lagrange-y, inf);
   
   n = points;
   x_lag_chebyshev = 1 - cos(pi*(0:n-1)/(n-1));
   y_lag_chebyshev = exp(x_lag_chebyshev);
   
   tic
   result_lag_cheby = polyinterp(x_lag_chebyshev,y_lag_chebyshev, x);
   lag_time_cheby(k) = toc;
   
   norm_2_lag_cheby(k) = norm(result_lag_cheby - y,2);
   norm_inf_lag_cheby(k) = norm(result_lag_cheby - y,inf);
   
   k = k+1;

end 

norm_2_bary = zeros(1,length(all_test_cases));
norm_inf_bary = zeros(1,length(all_test_cases));

norm_2_bary_cheby = zeros(1,length(all_test_cases));
norm_inf_bary_cheby = zeros(1,length(all_test_cases));


bary_time = zeros(1,length(all_test_cases));
bary_time_cheby = zeros(1,length(all_test_cases));


k = 1;
for points = all_test_cases
   x_bary = linspace(0,2,points);
   y_bary = exp(x_bary);
   
   tic
   result_baryrange = barylag([x_bary;y_bary]', x);
   bary_time(k) = toc;
   
   norm_2_bary(k) = norm(result_baryrange-y,2);
   norm_inf_bary(k) = norm(result_baryrange-y, inf);
   
   n = points;
   x_bary_chebyshev = 1 - cos(pi*(0:n-1)/(n-1));
   y_bary_chebyshev = exp(x_bary_chebyshev);
   
   tic
   result_bary_cheby = barylag([x_bary_chebyshev;y_bary_chebyshev]', x);
   bary_time_cheby = toc;
   
   norm_2_bary_cheby(k) = norm(result_bary_cheby - y,2);
   norm_inf_bary_cheby(k) = norm(result_bary_cheby - y,inf);
   
   k = k+1;

end 

%% 4 PCHIP

norm_2_pchip = zeros(1,length(all_test_cases));
norm_inf_pchip = zeros(1,length(all_test_cases));

norm_2_pchip_cheby = zeros(1,length(all_test_cases));
norm_inf_pchip_cheby = zeros(1,length(all_test_cases));


pchip_time = zeros(1,length(all_test_cases));
pchip_time_cheby = zeros(1,length(all_test_cases));

k = 1;
for points = all_test_cases
   x_pchip = linspace(0,2,points);
   y_pchip = exp(x_pchip);
   
   tic
   result_pchip = pchip(x_pchip,y_pchip, x);
   pchip_time(k) = toc;
   
   norm_2_pchip(k) = norm(result_pchip-y,2);
   norm_inf_pchip(k) = norm(result_pchip-y, inf);
   
   n = points;
   x_pchip_chebyshev = 1 - cos(pi*(0:n-1)/(n-1));
   y_pchip_chebyshev = exp(x_pchip_chebyshev);
   
   tic 
   result_pchip_cheby = pchip(x_pchip_chebyshev, y_pchip_chebyshev, x);
   pchip_time_cheby = toc;
   
   norm_2_pchip_cheby(k) = norm(result_pchip_cheby - y,2);
   norm_inf_pchip_cheby(k) = norm(result_pchip_cheby - y,inf);
   
   k = k+1;

end 