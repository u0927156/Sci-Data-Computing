function y = apply_operator_2d(v)
%
% applies discrete 2nd order approx to -lap(v)  assumes M=N
%
M = sqrt(length(v))-1; 
h = 1/M; hsq = h*h;
vm = reshape(v,M+1,M+1);
fac = 3.0;

ym = zeros(M+1,M+1);
ym(:,1) = 1/fac; ym(1,:) = 1/fac; 
for jj= 1:M+1
ym(M+1,jj) = exp(jj/(M+1))/fac; 
ym(jj,M+1) = exp(jj/(M+1))/fac;
end

ym(2:M,2:M) = ( (-vm(1:M-1,2:M)+2*vm(2:M,2:M)-vm(3:M+1,2:M)) ...
             +  (-vm(2:M,1:M-1)+2*vm(2:M,2:M)-vm(2:M,3:M+1)) )/(hsq);

% put ym in y 

y = ym(:);
