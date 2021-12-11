%% 1. Newton's method
clear; clc
f =@(x) [x(1).^3 + x(2) - 1; x(2).^3 - x(1) + 1];

df = @(x)[3*x(1).^2, 1;
    -1, 3.*x(2).^2];


% stop criterion
epsilon = 10^-10;

% Set up starting variables 
starting_x = [ 1.0, 0.1
    2.0 -.01
    100, 100
    0 0 
    0 100];

display('x1_0 & x2_0 & # iterations & x1 Solution & x2 Solution')
for n = 1:size(starting_x,1)    

    x1 = starting_x(n, 1);
    x2 = starting_x(n, 2);
    
    [x1_solution, x2_solution, num_iters] = ...
        newtons_method_2d(f, df, epsilon, x1, x2, 1000);
    
    display([num2str(x1) ' & ' num2str(x2) ' & ' num2str(num_iters)...
        ' & ' num2str(x1_solution) ' & ' num2str(x2_solution)])
    
end

%% Harder function
clc
f =@(x) [x(1).^2 + x(2).^2 - 2; exp(x(1)-1) + x(1)^2 - 3];

df = @(x)[2*x(1), 2*x(2);
    exp(x(1)), -2*x(2)];

% Set up starting variables 
starting_x = [ 1.29 0.58
    1.1, 1.1
    2.0 0.5
    3.0 5.0
    -0.7 1.14];

display('Normal Case')
display('Initial x1 & Initial x2 & num iterations & x1 Solution & x2 Solution')
for n = 1:size(starting_x,1)    

    x1 = starting_x(n, 1);
    x2 = starting_x(n, 2);
    
    [x1_solution, x2_solution, num_iters] = ...
        newtons_method_2d(f, df, epsilon, x1, x2, 1000);
    
    
    display([num2str(x1) ' & ' num2str(x2) ' & ' num2str(num_iters)...
        ' & ' num2str(x1_solution) ' & ' num2str(x2_solution)])
    
end
 %% 
 clc
display('Swapped to -2 in second equation')
display('Initial x1 & Initial x2 & num iterations & x1 Solution & x2 Solution')
f =@(x) [x(1).^2 + x(2).^2 - 2; exp(x(1)-1) + x(1)^2 - 4];
for n = 1:size(starting_x,1)    

    x1 = starting_x(n, 1);
    x2 = starting_x(n, 2);
    
    [x1_solution, x2_solution, num_iters] = ...
        newtons_method_2d(f, df, epsilon, x1, x2, 1000);
   
    display([num2str(x1) ' & ' num2str(x2) ' & ' num2str(num_iters)...
        ' & ' num2str(x1_solution) ' & ' num2str(x2_solution)])
    
end

%% Loran
clc
f = @(x) [x(1).^2/186^2 - x(2).^2/(300^2-186^2) - 1;
    (x(2)-500).^2/279^2 - (x(1)-300).^2/(500^2-279^2) - 1];


df = @(x)[ 2*x(1)/186^2 , - 2 * x(2)/(300^2-186^2);...
    (2 * (-300 + x(1) ) )./172159, (2 * (-500 + x(2) ) )/77841];

starting_x = [ 400:600; 400:600]';
display('x1_0 & x2_0 & # iterations & x1 Solution & x2 Solution')
for n = 1:size(starting_x,1)    

    x1 = starting_x(n, 1);
    x2 = starting_x(n, 2);
    
    [x1_solution, x2_solution, num_iters] = ...
        newtons_method_2d(f, df, epsilon, x1, x2, 25);
   
    display([num2str(x1) ' & ' num2str(x2) ' & ' num2str(num_iters)...
        ' & ' num2str(x1_solution) ' & ' num2str(x2_solution) '\\'])
    
end


%% 4 8.21 in Holmes
clc; clear
x = -2:.1:2;
y = -2:.1:2;

[X Y] = meshgrid(x, y);

Fxy = 1/10 .* (X+Y).^4 + (X-1).^2 + 4 .*Y.^2;

f = @(x)(1/10 .* (x(1)+x(2)).^4 + (x(1)-1).^2 + 4 .*x(2).^2);

surf(X, Y, Fxy)

% Steepest Descent Method Cauchy Method 
% 4(x-y)^3+4x-1 =0  ,  4(x-y)^3+2y+2
n=0;            %initialize iteration counter
eps=1;          %initialize error
a=0.09;         %set iteration parameter
x=[1;1];        %set starting value
%Computation loop
display('n & epsilon & x & y \\')
while eps>1e-10&n<100
    gradf=[2 * (-1 + x(1)) + 2/5 * (x(1) + x(2))^3;
    8 * x(2) + 2/5 * (x(1) + x(2))^3];  %gradf(x)
    eps=abs(gradf(1))+abs(gradf(2));                             %error
    y=x-a*gradf;                                                 %iterate
    x=y;                                                         %update x
    n=n+1;     %counter+1
    fprintf(' %i & %9.3e & %9.3e & %9.3e\\\\ \n',n,eps,x) 
end
n,x,eps,        %display end values

%% Rosenbock Problem
clc

a = 1; b = 100;

R = (a-X).^2 + b * (Y - X.^2).^2;

surf(X,Y,R)
title('Rosenbrock Function')
xlabel('X'); ylabel('Y')

% Steepest Descent Method Cauchy Method 

eps=1;          %initialize error
a=0.001;         %set iteration parameter



starting_x = [ 0.5 0.5
    0.1 0.1
    -.5 .5];


%Computation loop
display('Rosenbrock Problem')
    display('$x_0$ & $y_0$ & n & eps & x & y')
for m = 1:size(starting_x,1) 
    x0 = starting_x(m,:)';
    x = starting_x(m,:)';   %set starting value
    n=0;         %initialize iteration counter
    eps=1;
    
    
    while eps>1e-8 & n<100000
        gradf=[400*x(1)^3-400*x(1)*x(2)+2*x(1)-2;
        200*(x(2)-x(1)^2)]; %gradf(x)
        eps=abs(gradf(1))+abs(gradf(2));                             %error
        y=x-a*gradf;                                                 %iterate
        x=y;                                                         %update x
        n=n+1;     %counter+1
        %fprintf(' n= %i & eps = %9.3e & x= %9.3e  & y=%9.3e \n',n,eps,x) 
    end
    fprintf('%9.3f & %9.3f & %i & %9.3e & %9.3f & %9.3f\\\\ \n',x0, n,eps,x) 
    %n,x,eps,        %display end values
end

%% Armijo Algorithm
clc 

f = @(x)((a-x(1)).^2 + b * (x(2) - x(1).^2).^2);

    
gamma = 1/1000;
tau = 2/4;
a = 1;

v = [0.5 0.5];
n = 1;
eps = 1;

starting_x = [ 0.5 0.5
    0.1 0.1
    -.5 .5];


%Computation loop

for m = 1:size(starting_x,1) 
    n = 0;
    v = starting_x(m,:) ;
    eps = 1;
    while eps>1e-8 & n<1000000
        a = 1;
        gradf=[400*v(1)^3-400*v(1)*v(2)+2*v(1)-2;
        200*(v(2)-v(1)^2)];  %gradf(x)
        eps=abs(gradf(1))+abs(gradf(2));
        g = gradf;
        d = -g;

        while ~(f(v + a * d') <= f(v) + a * gamma * dot(g, d)) & a ~= 0
            a  = tau * a;

        end

        v = v + (a * d)';

        n = n+1;
    end
    n
    eps
    gamma;
end

%% Adam
clc

theta = [0.5 0.5]'

alpha = 1;
beta1 = 0.9;
beta2 = 0.999;
epsilon = 10^-8;

m = 0;
v = 0;
t = 0;

n = 0;
    
eps = 1;
while eps>1e-8 & n<1000000

    gradf=[400*theta(1)^3-400*theta(1)*theta(2)+2*theta(1)-2;
        200*(theta(2)-theta(1)^2)];

    g = gradf;

    m = beta1 * m + (1-beta1) * g;
    v = beta2 * v + (1 - beta2) * g.^2;

    m_hat = m / (1 - beta1);
    v_hat = v / (1 - beta2);

    theta = theta - alpha * m_hat ./(sqrt(v_hat) + epsilon);
    n = n+1;
end

theta
n
%% Function Definitions
[x1_solution, x2_solution, num_iters] = ...
    newtons_method_2d(f, df, epsilon, x1, x2)
function [x1, x2, num_iters] = newtons_method_2d(f, df, epsilon, x10, x20, max_iters)
x1 = x10;
x2 = x20;
num_iters = 0;
while norm(f([x1,x2])) > epsilon & num_iters < max_iters
    new_x = [x1; x2] - df([x1, x2])^-1 * f([x1,x2]);
    
    x1 = new_x(1);
    x2 = new_x(2);
    
    
    num_iters = num_iters + 1;
end

end