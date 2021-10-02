%% Problem 1
clear;
% Enter parameters
DeltaX = 1; H1 = 8; Hr =4;
as = [1.0 1.0e-5 1.0e-10, 1.0e-15]; % Make range of a 
n =161; % Select size of matrix
k = 1; % counter variable

conds = zeros(1,length(as)); % Preallocate matrix for condition numbers

for a = as
    % create the matrix with current a
    A = WaterMatrix(DeltaX, a, n);
    
    % Calculate condition number
    conds(k) = cond(A);
    
    % set up b
    b = zeros(n,1);
    b(1) = -H1; b(end) = -a*Hr;

    % select number of iterations for each method
    iterations = 400;

    % Calculate solution using Jacobi Method
    [x_jacobi, error_jacobi] = Jacobi_SP(A,b,iterations);

    % Calculate solution using Gauss-Seidel Method
    [x_gs, error_gs] = Gauss_Seidel_SP(A,b,iterations);

    % Calculate solution using Successive Over Relaxation
    omega = 1.95;
    [x_SOR, error_SOR] = SOR_SP(A,b,omega, iterations);

    % Plot results for each a
    figure(100+k)
    hold off
    semilogy(error_jacobi)
    hold on
    semilogy(error_gs, '--')
    semilogy(error_SOR, '-.k')

    xlabel('Iterations'); ylabel('Norm_2 Error')
    legend('Jacobi', 'Gauss-Seidel', 'SOR')
    title(['Error of Solutions of Water Flow Matrix a=' num2str(a)])
    
    saveas(gcf, ['10', num2str(k), '.png'])
    k = k + 1;
end

%% Plot Condition Number Trend
figure(105)
loglog(as, conds, '--ok')
set ( gca, 'xdir', 'reverse' ) % it makes more sense to look at how it changes 
% as a gets smaller. 
xlabel('a'); ylabel('Condition Number')
title('Condition Number of Flow of Water Matrix')
saveas(gcf, '105.png')
%% Problem 2 i

% times and errors recorded from 2DMGcode/V2d_driver
times = [.0566    0.0246    0.1794    0.7131    2.6199...
    10.9340   52.6200  626.8058];
errors = [0.0003    0.0001    0.0001    0.0004    0.0009...
    0.0018    0.0036    0.0072];
trial_sizes = 2.^[7:14];

% plot results
figure(201)
semilogy(log2(trial_sizes), times)
xlabel('Log_2 of Trial Size'); ylabel('Time')
title('Time of Multigrid Method over Trial Size')
saveas(gcf, '201.png')

figure(202)
semilogy(log2(trial_sizes), errors)
xlabel('Log_2 of Trial Size'); ylabel('norm_\infty error')
title('Norm Inf Error of Multigrid Method over Trial Size')
saveas(gcf, '202.png')

%% Problem 2 ii

% Information taken from Laplace2d_Modified
grid_size_log = [6 7 8 9 10];
    
Jacobi_Grid_Error = [0.00072801 0.00026193 0.00081935 0.0017643, 0.003591];
Jacobi_Grid_Time = [0.903295 0.995299 5.730664  63.056762, 774.645135];
RedBlack_Grid_Error = [0.031298 0.015846 0.0078258 0.0039071 0.0060559];
RedBlack_Grid_Time = [ 0.110093 0.103791  0.128461 0.186855, 0.441689];

% Plot results of experiments
figure(211)
hold off
semilogy(grid_size_log, Jacobi_Grid_Time)
hold on
semilogy(grid_size_log, RedBlack_Grid_Time)
xlabel('Log_2 of Trial Size'); ylabel('Time (sec)')

title('Time of Jacobi and Red Black over Trial Size')
saveas(gcf, '211.png')

figure(212)
hold off
semilogy(grid_size_log, Jacobi_Grid_Error)
hold on
semilogy(grid_size_log, RedBlack_Grid_Error)
xlabel('Log_2 of Trial Size'); ylabel('norm_\infty error')

title('Norm_\infty Error of Jacobi and Red Black over Trial Size')
saveas(gcf, '212.png')
