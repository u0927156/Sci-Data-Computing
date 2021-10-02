function y = relax_2d(v,f,nu)
%
% 2D weighted Jacobi, w = 4/5 : assumes M=N
%
M = sqrt(length(v))-1; h = 1/M; hsq = h*h;
vm = reshape(v,M+1,M+1);  fm = reshape(f,M+1,M+1);

ym = vm;

% relax and set homogeneous Dirichlet BC

w = 4/5;
for k=1:nu
    ym(2:M,2:M) = (1-w)*ym(2:M,2:M) + ...
                  (w/4)*(ym(1:M-1,2:M)+ym(3:M+1,2:M)+ ...
                         ym(2:M,1:M-1)+ym(2:M,3:M+1)+ ...
                         hsq*fm(2:M,2:M));
%    ym(:,1) = 0; ym(1,:) = 0; ym(M+1,:) = 0; ym(:,M+1) = 0;
end

% put ym in y
y = ym(:);
