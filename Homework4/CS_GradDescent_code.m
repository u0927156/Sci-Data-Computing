% Steepest Descent Method Cauchy Method 
% 4(x-y)^3+4x-1 =0  ,  4(x-y)^3+2y+2=0
n=0;            %initialize iteration counter
eps=1;          %initialize error
a=0.09;         %set iteration parameter
x=[1;1];        %set starting value
%Computation loop
while eps>1e-10&n<100
    gradf=[4*(x(1)-x(2))^3+4*x(1)-1;-4*(x(1)-x(2))^3+2*x(2)+2];  %gradf(x)
    eps=abs(gradf(1))+abs(gradf(2));                             %error
    y=x-a*gradf;                                                 %iterate
    x=y;                                                         %update x
    n=n+1;     %counter+1
    fprintf(' n= %i eps = %9.3e x(1)= %9.3e x(2)=%9.3e \n',n,eps,x) 
end
n,x,eps,        %display end values

