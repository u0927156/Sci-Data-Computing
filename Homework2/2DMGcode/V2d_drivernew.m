%
%  2D V-cycle for -(u_xx + u_yy) = f 
%                   u_{\Gamma} = 0   (homogeneous BC)
%
%  setup problem 
%
M = 128; hx = 1/M;  x = (0:hx:1)';
N = 128; hy = 1/N;  y = (0:hy:1)';
u = zeros(M+1,N+1);
f = u;


for i = 1:M+1
    for j = 1:M+1
        u(i,j) = exp(x(i)*y(j))/3.0;
        f(i,j) = u(i,j)*(x(i)^2+y(j)^2);
    end
end
%u(1,:) = 0;  u(M+1,:) = 0;  u(:,1) = 0;  u(:,M+1) = 0;

% normalize u and set f

%fac = 8.0*exp(2); 
%f = f/fac; 
%u = u/fac;

% setup call to V-cycle

v_in =  0*u;
nu_down = 5;
nu_up   = 5;
tic
v_out = Vcycle_2d(v_in,f,nu_down,nu_up);
toc
v_out = reshape(v_out,M+1,M+1);
norm(u-v_out,'inf')



