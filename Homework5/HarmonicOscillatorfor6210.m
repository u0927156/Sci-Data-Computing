% Example: Euler's method amd symplectic euler for harmonic oscillator 
% NEED TO REWRITE AS HARMONIC OSCILLATOR IN TWO AND HREE EQUATION FORM
clear all
clc

for nruns = 1:1
%********************************************************
%
if(nruns == 1)dt = 1e-3;end
%
omega = 100.0;
%if(nruns == 1)dt = 0.004;end
ts = 20.0;
npoints = int64(ts/dt);% standard  euler
ym1 = zeros(npoints,1); %x 
y0 = zeros(npoints,1);  %v
 
% standard symplectic euler A
y1 = zeros(npoints,1); % x
y2 = zeros(npoints,1); % v 

% Stormer Verlet 
x_sv = zeros(npoints,1); % x
v_sv = zeros(npoints,1); % v 
H_sv = zeros(npoints,1);

% Reverse Euler
x_re = zeros(npoints,1); % x
v_re = zeros(npoints,1); % v 
H_re = zeros(npoints,1);

% Trapezoidal 
x_tr = zeros(npoints,1); % x
v_tr = zeros(npoints,1); % v 
H_tr = zeros(npoints,1); 


xtrue = zeros(npoints,1);
vtrue = zeros(npoints,1);

Htrue = zeros(npoints,1);

t = zeros(npoints,1);

X0 = 1.0/sqrt(2.0);
%X0 = 1.0;
V0 = sqrt( 1.0 - X0^2);
 
xtrue(1) = X0;
vtrue(1) = V0;
Htrue(1) = 1/2 * (vtrue(1)^2 + omega^2 * xtrue(1)^2);

ym1(1) = X0; % the initial condition
y0(1) =  V0;

y1(1) = ym1(1); % the initial condition
y2(1) = y0(1);

x_sv(1) = X0; % Stormer Verlet Initial Conditions
v_sv(1) = V0;
H_sv(1) = 1/2*(v_sv(1)^2 + omega^2 * x_sv(1)^2);


x_re(1) = X0; % Reverse Euler Initial Conditions
v_re(1) = V0;
A = [1 -dt; dt*omega^2 1]; % Matrix to be inverted for implicit Euler's

H_re(1) = 1/2*(v_re(1)^2 + omega^2 * x_re(1)^2);

Ainv = inv(A);

x_tr(1) = X0; % Trapezoidal Initial Conditions
v_tr(1) = V0;
H_tr(1) = 1/2*(v_tr(1)^2 + omega^2 * x_tr(1)^2);

Binv = inv([1 -1/2*dt; 1/2*dt*omega^2 1]);
C = [1 1/2*dt; -1/2*dt*omega^2 1];




t(1) = 0.0;
%
for step=1:npoints-1 % loop over the timesteps
    % standard Euler
    % ym1 is x  and y0 is y
    ym1(step+1) = ym1(step) + dt*y0(step);
    y0(step+1) = y0(step) - dt*omega^2 *ym1(step);
            
    %symplectic Euler Euler A
    % y1 is x  and y2 is y
    y1(step+1) = y1(step) + dt*y2(step);
    y2(step+1) = y2(step) - dt* omega^2 *y1(step+1);

    % Stormer-Verlet Method
    % x_sv is x  and v_sv is y
    v_half = v_sv(step) - dt/2 *  omega^2 *  x_sv(step);
    x_sv(step+1) = x_sv(step) + dt * v_half;
    v_sv(step+1) = v_half - dt/2 * (omega^2 * x_sv(step+1));
    H_sv(step+1) = 1/2*(v_sv(step+1)^2 + omega^2 * x_sv(step+1)^2);
    
    % Reverse Euler Method
    temp = Ainv * [x_re(step); v_re(step)];
    x_re(step+1) = temp(1);
    v_re(step+1) = temp(2);
    H_re(step+1) = 1/2*(v_re(step+1)^2 + omega^2 * x_re(step+1)^2);
    
    % Trapezoidal
    temp = Binv * C * [x_tr(step); v_tr(step)];
    x_tr(step+1) = temp(1);
    v_tr(step+1) = temp(2);
    H_tr(step+1) = 1/2*(v_tr(step+1)^2 + omega^2 * x_tr(step+1)^2);
    
    
      t(step+1) = t(step) + dt;
% Exact 
     a = cos(omega*t(step+1)); b = sin(omega*t(step+1));
     xtrue(step+1) = a *xtrue(1) + b/omega*vtrue(1);
     vtrue(step+1) = -b*omega*xtrue(1) + a*vtrue(1);
     Htrue(step+1) = 1/2 * (vtrue(step+1)^2 + omega^2 * xtrue(step+1)^2);
     
     ii = step+1;

end
%final values 
fprintf('xt=%5.2e x1 %5.2e x2=%5.2e ',xtrue(ii),ym1(ii),y1(ii))
fprintf('yt=%5.2e v1 %5.2e v2=%5.2e ',vtrue(ii), y0(ii),y2(ii))
fprintf('\n')
%*************************************************************************************
  figure(4)
  subplot(1,2,1)
  plot(ym1,y0,xtrue,vtrue)  %plots x and v 
  legend(' Euler ', 'true')
  xlabel('x ') % x-axis label
  ylabel('v ') % y-axis label
  
    subplot(1,2,2)
  plot(y1,y2,xtrue,vtrue)  %plots x and v 
%  legend(' omega = 2.0 .004','y2 dt=0.004','Exacty1','Exacty2')
  legend(' Symplectic Euler A ', 'true')
  xlabel('x ') % x-axis label
  ylabel('v ') % y-axis label

   figure(5)
   
   plot(x_sv, v_sv, '.-')
   hold on
   plot(xtrue, vtrue, '--')
   hold off
    legend('Stormer Verlet', 'true')
  xlabel('x ') % x-axis label
  ylabel('v ') % y-axis label
  
  figure(6)
   plot(x_re, v_re)
   hold on
   plot(xtrue, vtrue, '--')
   hold off
    legend('Reverse Euler', 'true')
  xlabel('x ') % x-axis label
  ylabel('v ') % y-axis label
  
  figure(7)
   plot(x_tr, v_tr, '.-')
   hold on
   plot(xtrue, vtrue, '-.')
   hold off
    legend('Trapezoidal', 'true')
  xlabel('x ') % x-axis label
  ylabel('v ') % y-axis label
end

display(norm(Htrue - H_sv, 'inf'))
display(norm(Htrue - H_re, 'inf'))
display(norm(Htrue - H_tr, 'inf'))

