function y = apply_operator_2d(v)
%
% applies discrete 2nd order approx to -lap(v)  assumes M=N
%
M = sqrt(length(v))-1; 
h = 1/M; hsq = h*h;
vm = reshape(v,M+1,M+1);

ym = zeros(M+1,M+1);
ym(2:M,2:M) = ( (-vm(1:M-1,2:M)+2*vm(2:M,2:M)-vm(3:M+1,2:M)) ...
             +  (-vm(2:M,1:M-1)+2*vm(2:M,2:M)-vm(2:M,3:M+1)) )/(hsq);
ym(:,1) = 0; ym(1,:) = 0; ym(M+1,:) = 0; ym(:,M+1) = 0;

% put ym in y 

y = ym(:);
