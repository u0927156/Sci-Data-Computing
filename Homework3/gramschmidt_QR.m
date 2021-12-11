function [Q,R] =gramschmidt_QR(A)
Q=A;
n = size(Q,2);
R = zeros(n,n);
for j = 1:n
    q = A(:,j);
    for k = 1:j-1
        qk= Q(:,k);
        R(k,j)= (qk'*q);
        q = q -qk*R(k,j);
    end
    R(j,j)= norm(q,2);
    Q(:,j)=q/R(j,j);
end
