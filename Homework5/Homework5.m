%% Problem 1, Simpson's rule, polynomials
% Finds the area underneath the function x^p using Simpson's rule

clear; clc


ps = [2 3 4 5 6 8]; % The values of p
Ns = [17 33 65 129 257 513]; % The number of intervals to calculate with

a = 0; b = 1; % The range to calculate on

% Store values for display later.
computed_values = zeros(length(ps), length(Ns));


for m = 1:length(ps)
    p = ps(m); % Get current p
    F = @(x) (x.^p); % Make the function
    for n = 1:length(Ns)
        N = Ns(n); % Get current N


        % Calculate area.
        computed_values(m,n) = simpsons_rule(F, a, b, N);

    end
end

computed_values

%% Problem 1, Simpson's Rule, Sinusoidal Function
clc
F = @(x) (1 + sin(x) .* cos(2.*x./3) .* sin(4.*x));
a = 0; b = 2*pi;

Ns = [ 17 33 65 129 257 513];

sinusoidal_function_values = zeros(1, length(Ns));
for m = 1:length(Ns)
   N = Ns(m);
   
   sinusoidal_function_values(m) = simpsons_rule(F, a,b, N);
    
end
sinusoidal_function_values
%% Problem 2, QuadTX 
clc; 

F = @(x) (cos(x.^3).^200);
a = 0; b = 3;

tols = 10.^-(7:14);

counts = zeros(1, length(tols));
areas = zeros(1, length(tols));
times = zeros(1, length(tols));

for m = 1:length(tols)
    tol = tols(m);
    tic
    [Q,fcount] = quadtx(F,a,b, tol);
    times(m) = toc;
    areas(m) = Q;
    counts(m) = fcount;
end

areas
counts
times

temp = sprintf(' %.0s&', tols);
temp(end) = [];             %get rid of trailing comma
fprintf(temp)
fprintf('\n')

temp = sprintf(' %d&', counts);
temp(end) = [];             %get rid of trailing comma
fprintf(temp)
fprintf('\n')

temp = sprintf(' %.2f ms &', times.*1000);
temp(end) = [];             %get rid of trailing comma
fprintf(temp)
fprintf('\n')

%% Problem 3
clc
Ns = 2.^(11:17)+1;

simpsons_tough = zeros(1, length(Ns));
simpsons_times = zeros(1, length(Ns));
for m = 1:length(Ns)
    N = Ns(m);
    tic
    simpsons_tough(m) = simpsons_rule(F, a,b, N);
    simpsons_times(m) = toc;
end

simpsons_tough

tic
[Q,fcount] = quadtx(F,a,b, tol);
quadtx_time = toc;

%% Problem 4, Cooling 
clear; clc
minutes = 5;
total_time= minutes*60;
t = 0:.01:(total_time);
r = 0.025; % s^-1 
Ts = 19; % Degrees Celcius

T_analytical = 65.*exp(-r.*t) + Ts;
T = @(t) (65.*exp(-r.*t) + Ts);

figure(1)
plot(t, T_analytical, '--')
hold on;
hs = [30 15 10 5 1 .5 .25];

first_step_errors = ones(1, length(hs));
for m = 1:length(hs)
    
    h = hs(m);
    
    curr_t = 0;

    dT = @(Tc) (-r .* (Tc - Ts));
    T0 = 84;

    curr_temp = T0;
    temps = [];
    times = [];
    while curr_t < total_time
        slope = dT(curr_temp);

        next_temp = curr_temp + h * slope;
        next_t = curr_t + h;

        times = [times curr_t];
        temps = [temps curr_temp];

        curr_t = next_t;
        curr_temp = next_temp;

    end

    plot(times, temps, '-')
    
    after_first_step = times(2);
    T_analytical_after_first_step = T(after_first_step);
    T_estimated_after_first_step = temps(2);
    
    first_step_errors(m) = T_analytical_after_first_step - T_estimated_after_first_step;
    
end
legend('Analytical', '30', '15', '10', '5', '1', '0.5', '0.25')
hold off
title('Forward Euler Solution of Coffee Cup Problem')
xlabel('Time (s)')
ylabel('Temperature (C^\circ)')




temp = sprintf(' %.1d&', hs);
temp(end) = [];             %get rid of trailing comma
fprintf(temp)
fprintf('\n')

temp = sprintf(' %.3d&', first_step_errors);
temp(end) = [];             %get rid of trailing comma
fprintf(temp)
fprintf('\n')
%% Problem 5-6, ODE23
clc;
f = @(t,y)(-r .* (y - Ts));


curr_t = 0;

y0 = 84;

figure(2)
plot(t, T_analytical, '--')
hold on;
hs = [30 15 10 5 1 .5 .25];

predicted_errors_ode23 = zeros(1, length(hs));
actual_errors_ode23 = zeros(1, length(hs));


for  m = 1:length(hs)
    h = hs(m);
    [times, ys, errors] = MyODE23(f, 0, 300, y0, h);

    plot(times, ys, '-.')
    
    predicted_errors_ode23(m) = errors(2);
    actual_errors_ode23(m) = T(times(2)) - ys(2) ;
end
legend('Analytical', '30', '15', '10', '5', '1', '0.5', '0.25')
hold off
title('ODE 23 Solution of Coffee Cup Problem')
xlabel('Time (s)')
ylabel('Temperature (C^\circ)')

temp = sprintf(' %.3d&', predicted_errors_ode23);
temp(end) = [];             %get rid of trailing comma
fprintf(temp)
fprintf('\n')

temp = sprintf(' %.3d&', actual_errors_ode23);
temp(end) = [];             %get rid of trailing comma
fprintf(temp)
fprintf('\n')
%% Problem 7, Bigger r

r = 0.6;
f = @(t,y)(-r .* (y - Ts));


curr_t = 0;

y0 = 84;

figure(3)
subplot(2,1,1)
plot(t, T_analytical, '--')
hold on;
hs = [30 15 10 5 1 .5 .25];

subplot(2,1,2)
hold on;
for  h = hs
    
    [times, ys, errors] = MyODE23(f, 0, 300, y0, h);
    subplot(2,1,1)
    plot(times, ys, '-.')
    
    subplot(2,1,2)
    plot(times, errors, '-.')
end
subplot(2,1,1)
legend('Analytical', '30', '15', '10', '5', '1', '0.5', '0.25', 'location', 'Northwest')
hold off
title('ODE 23 Solution of Coffee Cup Problem, Big r')
xlabel('Time (s)')
ylabel('Temperature (C^\circ)')

subplot(2,1,2)
legend('30', '15', '10', '5', '1', '0.5', '0.25', 'location', 'Northwest')
title('ODE 23 Error of Coffee Cup Problem, Big r')
xlabel('Time (s)')
ylabel('Error')

hold off