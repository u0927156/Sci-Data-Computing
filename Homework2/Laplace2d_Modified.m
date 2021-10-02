% Solving the 2-D Laplace's equation by the Finite Difference
...Method 
% Numerical scheme used is a second order central difference in space
...(5-point difference)
clear;


%% Updated parameters, Used for both methods
%Specifying parameters
pow = 10;
nx=2^pow;                           %Number of steps in space(x)
ny=2^pow;                           %Number of steps in space(y)       
niter=1000;                     %Maximum Number of iterations 
dx=1/(nx-1);                     %Width of space step(x)
dy=1/(ny-1);                     %Width of space step(y)
x=0:dx:1;                        %Range of x(0,1) and specifying the grid points
y=0:dy:1;                        %Range of y(0,) and specifying the grid points
TOL=1e-18;                       % Tolerance to decide when to break

%% Calculate Actual Solution for Calculating Error

[X, Y] = meshgrid(x,y);
actual_solution = (sin(88*pi.*X)+sin(72*pi.*Y))./ ((88*pi)^2 + (72*pi)^2);

%% Set up initial conditions
%Initial Conditions
p=zeros(ny,nx);                  %Preallocating p, x_k+1
pn=zeros(ny,nx);                 %Preallocating pn, x_k
pnn = zeros(ny,nx);              %Preallocating pnn, which equivalent to x_k-1
%%
%Boundary conditions
p(:,1) = u_true(0,y); % y = y, x = 0 
p(:,nx)= u_true(x(end), y); % y = y, x = 1
p(1,:) = u_true(x, 0); % y=0, x = x                   
p(ny,:)= u_true(x,y(end));               

%% Modified Solution to solve the more difficult problem
%Explicit iterative scheme with C.D in space (5-point difference)
j=2:nx-1;
i=2:ny-1;
[J, I] = meshgrid(j,i); % need meshgrid of values that are changing for 
% solution
tic
for it=1:niter
    pnn = pn;
    pn=p;
    p(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1))))/...
        (2*(dx^2+dy^2))...
        + (sin(88*pi.*x(J)) + sin(72*pi.*y(I)))*((dx/4)^2);
    %Boundary conditions 
    p(:,1) = u_true(0,y); % y = y, x = 0 
    p(:,nx)= u_true(x(end), y); % y = y, x = 1
    p(1,:) = u_true(x, 0); % y=0, x = x                   
    p(ny,:)= u_true(x,y(end));  % y = 1, x=x 
    
    % Don't check tolerance with no previous information.
    if it == 1
        continue
        
    end
    
    % Calculate C
    c = norm(p - pn, 2) / norm(pn-pnn, 2);
    
    % Break if difference is below tolerance.
    if c/(1-c) * norm(p-pn,2) <= TOL
       break; 
    end
 end


%Plotting the solution
figure(301)
surf(x,y,p,'EdgeColor','none');       
shading interp
title({'Multigrid Problem, Solved with Jacobi Method';['{\itNumber of iterations} = ',num2str(it)]})
xlabel('Spatial co-ordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial co-ordinate (y)')
zlabel('Solution profile (P) \rightarrow')

display(['Jacobi Norm_inf Error =' num2str(norm(p-actual_solution,'inf'))]);
display(['Number of Iterations =' num2str(it)]);
toc
%% Red-Black Version

% Reset the variables
%Initial Conditions
p=zeros(ny,nx);                  %Preallocating p, x_k+1
pn=zeros(ny,nx);                 %Preallocating pn, x_k
pnn = zeros(ny,nx);              %Preallocating pnn, which equivalent to x_k-1

%Boundary conditions
p(:,1) = u_true(0,y); % y = y, x = 0 
p(:,nx)= u_true(x(end), y); % y = y, x = 1
p(1,:) = u_true(x, 0); % y=0, x = x                   
p(ny,:)= u_true(x,y(end));   



j=2:nx-1;
i=2:ny-1;
[J, I] = meshgrid(j,i); % need meshgrid of values that are changing for 
% solution

% Get indices for odd and even cases
i_even = i(1:2:end);
i_odd  = i(2:2:end);

j_even = j(1:2:end);
j_odd  = j(2:2:end);
    
% get meshgrid values for even and odd cases
[J_even, I_even] = meshgrid(j_even, i_even);
[J_odd, I_odd] = meshgrid(j_odd, i_odd);
tic
for it=1:niter
    % update previously stored info.
    pnn = pn;
    pn=p;
    
    % Red = Even
    p(i_even,j_even)=((dy^2*(pn(i_even+1,j_even)+pn(i_even-1,j_even)))+...
        (dx^2*(pn(i_even,j_even+1)+pn(i_even,j_even-1))))/...
        (4*(dx^2+dy^2))...
        + (sin(88*pi.*x(J_even)) + sin(72*pi.*y(I_even)))*((dx*2)^2);
    
    % Black = odd
    p(i_odd,j_odd)=((dy^2*(p(i_odd+1,j_odd)+p(i_odd-1,j_odd)))+...
        (dx^2*(p(i_odd,j_odd+1)+p(i_odd,j_odd-1))))/...
        (4*(dx^2+dy^2))...
        + (sin(88*pi.*x(J_odd)) + sin(72*pi.*y(I_odd)))*((dx*2)^2);
    
    %Boundary conditions 
    p(:,1) = u_true(0,y); % y = y, x = 0 
    p(:,nx)= u_true(x(end), y); % y = y, x = 1
    p(1,:) = u_true(x, 0); % y=0, x = x                   
    p(ny,:)= u_true(x,y(end));  % y = 1, x=x 
    
    % Skip tolerance check on first iteration because we have no previous
    % information
    if it == 1
        continue
        
    end
    
    % calculate c
    c = norm(p - pn, 'inf') / norm(pn-pnn, 'inf');
    
    if it > 100
        break
    end
    % Break if difference is below tolerance.
    if c/(1-c) * norm(p-pn,'inf') <= TOL
       break; 
    end
 end


%Plotting the solution
figure(302)
surf(x,y,p,'EdgeColor','none');       
shading interp
title({'Multigrid Problem, Solved with Red-Black Method';['{\itNumber of iterations} = ',num2str(it)]})
xlabel('Spatial co-ordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial co-ordinate (y)')
zlabel('Solution profile (P) \rightarrow')

display(['Red-Black Norm_inf Error =' num2str(norm(p-actual_solution,'inf'))]);
display(['Number of Iterations =' num2str(it)]);
toc
%% Define a function with the true value to make my life easier
function u = u_true(x,y)
    u = (sin(88*pi.*x)+sin(72*pi.*y))./ ((88*pi)^2 + (72*pi)^2);
end