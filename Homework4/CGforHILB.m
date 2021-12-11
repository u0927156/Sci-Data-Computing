
%% Conjugate Gradient example 
n = 4;
A = hilb(n)
b=ones(n,1)
x=zeros(n,1)
r = -A*x+b;
tol = 0.1e-2;
[x, niters] = cgsolve(A,b);
 resid = A*x-b;
 normRes= norm(resid,inf);
fprintf(' %i Resnrm  %8.2e  \n',niters,normRes)
fprintf(' x(1)= %8.2e x(2)= %8.2e x(3)= %8.2e x(4) = %8.2e\n',x(1),x(2),x(3),x(4))


