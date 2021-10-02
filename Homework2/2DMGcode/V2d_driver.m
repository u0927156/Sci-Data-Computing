%
%  2D V-cycle for -(u_xx + u_yy) = f 
%                   u_{\Gamma} = 0   (homogeneous BC)
%
%  setup problem 
%
clc
times = zeros(1, 8)
errors = zeros(1,8)
k = 1
for power = 7:14
    M = 2^power; hx = 1/M;  x = (0:hx:1)';
    N = 2^power; hy = 1/N;  y = (0:hy:1)';
    u = zeros(M+1,N+1);  


    for i = 1:M+1
        for j = 1:M+1
            u(i,j) = sin(88*pi*x(i))*sin(72*pi*y(j));
        end
    end
    u(1,:) = 0;  u(M+1,:) = 0;  u(:,1) = 0;  u(:,M+1) = 0;

    % normalize u and set f

    fac = ((88*pi)^2+(72*pi)^2);
    %fac = 1.0;
    f = u; 
    u = u/fac;

    % setup call to V-cycle

    v_in =  0*u;
    nu_down = 5;
    nu_up   = 5;
    tic
    v_out = Vcycle_2d(v_in,f,nu_down,nu_up);
    times(k) = toc
    v_out = reshape(v_out,M+1,M+1);
    errors(k) = norm(u-v_out,'inf')
    
    k = k + 1;
end


