%% Problem 1
clear; clc; close all
x = linspace(0,2, 1001); %% 1001 spaces between 0 and 2
y = exp(x); % Create the exponential function

figure(101) % plot it.
plot(x,y)
title('Simple Exponential Function')

%% Problem 2
clc



all_test_cases = [6 11 21 41 81 161 321 641]; % the number of points to use

k = 1; % counter for the norms and time

% Create vectors for storing the norm
norm_2_even = zeros(1,length(all_test_cases));
norm_inf_even = zeros(1,length(all_test_cases));

norm_2_cheby = zeros(1,length(all_test_cases));
norm_inf_cheby = zeros(1,length(all_test_cases));

% Create vectors for time.
vander_time = zeros(1,length(all_test_cases));
vander_time_cheby = zeros(1,length(all_test_cases));

for points = all_test_cases
   x_sampled = linspace(0,2,points); % Create evenly spaced points to test
   y_sampled = exp(x_sampled);
   
   tic % for timing
   V = vander(x_sampled); % make the Vandermonde matrix
   
   coeffs = V/y_sampled; % calculate the polynomial coefficients
    
   
   result = polyval(coeffs,x); % interpolate based on the coefficients
   vander_time(k) = toc; 
   
   % Get the error norms
   norm_2_even(k) = norm(result-y,2);
   norm_inf_even(k) = norm(result-y, inf);
   
   n = points;
   x_chebyshev = 1 - cos(pi*(0:n-1)/(n-1)); % get the chebyshev points
   y_chebyshev = exp(x_chebyshev);
   
   % Same as above
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

% Create vectors for norms and time
norm_2_lag = zeros(1,length(all_test_cases));
norm_inf_lag = zeros(1,length(all_test_cases));

norm_2_lag_cheby = zeros(1,length(all_test_cases));
norm_inf_lag_cheby = zeros(1,length(all_test_cases));

lag_time = zeros(1,length(all_test_cases));
lag_time_cheby = zeros(1,length(all_test_cases));

k = 1;
for points = all_test_cases
    
   % create the points and function
   x_lag = linspace(0,2,points);
   y_lag = exp(x_lag);
   
   % Interpolate
   tic
   result_lagrange = polyinterp(x_lag,y_lag, x); % use the polyinterp func
   lag_time(k) = toc;
   
   % Calculate norms
   norm_2_lag(k) = norm(result_lagrange-y,2);
   norm_inf_lag(k) = norm(result_lagrange-y, inf);
   
   % Do the same for Chebyshev points
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

% Create vectors for storing norms and time for barylag folder
norm_2_bary = zeros(1,length(all_test_cases));
norm_inf_bary = zeros(1,length(all_test_cases));

norm_2_bary_cheby = zeros(1,length(all_test_cases));
norm_inf_bary_cheby = zeros(1,length(all_test_cases));


bary_time = zeros(1,length(all_test_cases));
bary_time_cheby = zeros(1,length(all_test_cases));


k = 1;
for points = all_test_cases
   % Create the function
   x_bary = linspace(0,2,points);
   y_bary = exp(x_bary);
   
   tic
   result_baryrange = barylag([x_bary;y_bary]', x); 
   % This function needs input with two columns 
   bary_time(k) = toc;
   
   % Get the norms
   norm_2_bary(k) = norm(result_baryrange-y,2);
   norm_inf_bary(k) = norm(result_baryrange-y, inf);
   
   % Do the same for chebyshev points 
   n = points;
   x_bary_chebyshev = 1 - cos(pi*(0:n-1)/(n-1));
   y_bary_chebyshev = exp(x_bary_chebyshev);
   
   tic
   result_bary_cheby = barylag([x_bary_chebyshev;y_bary_chebyshev]', x);
   bary_time_cheby(k) = toc;
   
   norm_2_bary_cheby(k) = norm(result_bary_cheby - y,2);
   norm_inf_bary_cheby(k) = norm(result_bary_cheby - y,inf);
   
   k = k+1;

end 

%% 4 PCHIP

% Make the vectors for norms and time
norm_2_pchip = zeros(1,length(all_test_cases));
norm_inf_pchip = zeros(1,length(all_test_cases));

norm_2_pchip_cheby = zeros(1,length(all_test_cases));
norm_inf_pchip_cheby = zeros(1,length(all_test_cases));

pchip_time = zeros(1,length(all_test_cases));
pchip_time_cheby = zeros(1,length(all_test_cases));

k = 1;
for points = all_test_cases
    % Create the points and function
   x_pchip = linspace(0,2,points);
   y_pchip = exp(x_pchip);
   
   % Interpolate
   tic
   result_pchip = pchip(x_pchip,y_pchip, x);
   pchip_time(k) = toc;
   
   % Norms
   norm_2_pchip(k) = norm(result_pchip-y,2);
   norm_inf_pchip(k) = norm(result_pchip-y, inf);
   
   % Same for chebyshev points
   n = points;
   x_pchip_chebyshev = 1 - cos(pi*(0:n-1)/(n-1));
   y_pchip_chebyshev = exp(x_pchip_chebyshev);
   
   tic 
   result_pchip_cheby = pchip(x_pchip_chebyshev, y_pchip_chebyshev, x);
   pchip_time_cheby(k) = toc;
   
   norm_2_pchip_cheby(k) = norm(result_pchip_cheby - y,2);
   norm_inf_pchip_cheby(k) = norm(result_pchip_cheby - y,inf);
   
   k = k+1;

end 

%% Print tables to LaTeX. 

% Make a cell with all of the vectors to print to tables
all_the_tables_I_need_to_print = {
    {norm_2_even, norm_inf_even, 'Vandermonde Even-Spacing Interpolation'}
    {norm_2_cheby, norm_inf_cheby, 'Vandermonde Chebyshev-Spacing Interpolation'}
    {norm_2_lag, norm_inf_lag, 'LaGrange Even-Spacing Interpolation'}
    {norm_2_lag_cheby, norm_inf_lag_cheby, 'LaGrange Chebyshev-Spacing Interpolation'}
    {norm_2_bary, norm_inf_bary, 'Barylag Even-Spacing Interpolation'}
    {norm_2_bary_cheby, norm_inf_bary_cheby, 'Barylag Chebyshev-Spacing Interpolation'}
    {norm_2_pchip, norm_inf_pchip, 'Pchip Even-Spacing Interpolation'}
    {norm_2_pchip_cheby, norm_inf_pchip_cheby, 'Pchip Chebyshev-Spacing Interpolation'}
};

% Print all of the tables using a function below.
for n = 1:length(all_the_tables_I_need_to_print)
    PrintNorms(all_the_tables_I_need_to_print{n}{1},all_the_tables_I_need_to_print{n}{2}, all_the_tables_I_need_to_print{n}{3})
end
%

%% 5
% Make plots for time and errors

% Evenly spaced timing plot
figure(501)
semilogy(all_test_cases, [vander_time; lag_time; bary_time; pchip_time] ,'--*')
legend('Vandermonde', 'Lagrange', 'Barylag', 'PCHIP','location', 'northwest')
xlabel('Sampled Points'); ylabel('Time (s)'); title('Evenly Spaced Algorithm Time')
saveas(gcf, 'EvenTime.png')

% Chebyshev points timing plots
figure(502)
semilogy(all_test_cases, [vander_time_cheby; lag_time_cheby; bary_time_cheby; pchip_time_cheby] ,'--*')
legend('Vandermonde', 'Lagrange', 'Barylag', 'PCHIP', 'location', 'northwest')
xlabel('Sampled Points'); ylabel('Time (s)'); title('Chebyshev Points Algorithm Time')
saveas(gcf, 'ChebyTime.png')

% Error for evenly spaced points
figure(503)
semilogy(all_test_cases, [norm_inf_even; norm_inf_lag; norm_inf_bary; norm_inf_pchip] ,'--*')
legend('Vandermonde', 'Lagrange', 'Barylag', 'PCHIP', 'location', 'northwest')
xlabel('Sampled Points'); ylabel('Time (s)'); title('Infinity Norm Error for Even Spaced Points')
saveas(gcf, 'EvenError.png')

% Error for Chebyshev points
figure(504)
semilogy(all_test_cases, [norm_inf_cheby; norm_inf_lag_cheby; norm_inf_bary_cheby; norm_inf_pchip_cheby] ,'--*')
legend('Vandermonde', 'Lagrange', 'Barylag', 'PCHIP', 'location', 'northwest')
xlabel('Sampled Points'); ylabel('Time (s)'); title('Infinity Norm Error for Chebyshev Spaced Points')
saveas(gcf, 'ChebyError.png')

%% 6 Weather Interpolating


% load the data
load ud
load zd
load zp-1


% Vandermonde Interpolation
tic
V = vander(zd);
weather_c = V/ud';
result_vander = polyval(weather_c,zp_1); % Get the interpolated points at the physics points
V_back = vander(zp_1);
weather_back_c = V_back / result_vander'; % interpolate again
backresult_vander = polyval(weather_back_c, zd); % Get points back at dynamic points
times(1) = toc;



hold on;

% Lagrange interpolation
tic
result_lagrange = polyinterp(zd,ud, zp_1); % Interpolate at physics points
backresult_lagrange = polyinterp(zp_1, result_lagrange, zd);
% Interpolate back at dynamic points
times(2) = toc;



% Barycentric lagrange interpolation
tic
result_barylag = barylag([zd,ud], zp_1);
backresult_barylag = barylag([zp_1, result_barylag], zd);
times(3) = toc;



% Pchip interpolation. 
tic
result_pchip = pchip(zd,ud, zp_1);
backresult_pchip = polyinterp(zp_1, result_pchip, zd);
times(4) = toc;

% Plot the inteprolated results
figure(601)
plot(zd, ud)
hold on
plot(zd, backresult_vander)
plot(zd, backresult_lagrange)
plot(zd, backresult_barylag)
plot(zd, backresult_pchip)
hold off
legend('Original', 'Vander', 'Lagrange', 'Barylag', 'Pchip')
title('Weather Interpolation')

% Plot times to solve with each algorithm
figure(602)
hold on
plot(times, '*')
xticks([1 2 3 4])
xlim([.5 4.5])
xticklabels({'Vander', 'Lagrange', 'Barylag', 'Pchip'})
ylabel('Time (s)')
title('Algorithm Time for Weather Data')

% get the norms for each of the algoriths
norm_2 = zeros(1,4);
norm_inf = zeros(1,4);

norm_2(1) = norm(backresult_vander - ud,2);
norm_inf(1) = norm(backresult_vander - ud,inf);

norm_2(2) = norm(backresult_lagrange - ud,2);
norm_inf(2) = norm(backresult_lagrange - ud,inf);

norm_2(3) = norm(backresult_barylag - ud,2);
norm_inf(3) = norm(backresult_barylag - ud,inf);

norm_2(4) = norm(backresult_pchip - ud,2);
norm_inf(4) = norm(backresult_pchip - ud,inf);

% Print out the table with all of the errors.
row_format_2 = '2 & %4.3e & %4.3e &  %4.3e & %4.3e \\\\ \\hline \n';
row_format_inf = '$\\infty$ & %4.3e & %4.3e &  %4.3e & %4.3e \\\\ \\hline \n';
table_header = '\\hline norm & Vandermonde & Lagrange & Barylag & pchip \\\\ \\hline \\hline \n';
display('\begin{center}')
display('\begin{tabular}{| c | c | c | c | c |}')
fprintf(table_header)
fprintf(row_format_2, norm_2)
fprintf(row_format_inf, norm_inf)
display('\end{tabular}')
display('\end{center}')

%% Functions
function PrintNorms(norm_2, norm_inf, table_title)
   % Print the table for the norms with a given title.

    table_header = '\\hline norm & 6 & 11 & 21 & 41 & 81 & 161 & 321 & 641 \\\\ \\hline \\hline \n';
    
    row_format_2 = '2 & %4.2e & %4.2e &  %4.2e & %4.2e &  %4.2e & %4.2e &  %4.2e & %4.2e \\\\ \\hline \n';
    row_format_inf = '$\\infty$ & %4.2e & %4.2e &  %4.2e & %4.2e &  %4.2e & %4.2e &  %4.2e & %4.2e \\\\ \\hline \n'; %ideal_format = '\\hline %d & %d &  %d & %d &  %d & %d &  %d & %d \\\\ \\hline \\hline \n';
    %fprintf(ideal_format,ideal);
    
    fprintf(['\n' table_title '\n'])
    display('\begin{center}')
    display('\begin{tabular}{| c | c | c | c | c | c | c | c | c |}')
    fprintf(table_header)
    fprintf(row_format_2, norm_2)
    fprintf(row_format_inf, norm_inf)
    display('\end{tabular}')
    display('\end{center}')
    %fprintf(row_format,['\infty', norm_inf])
end
