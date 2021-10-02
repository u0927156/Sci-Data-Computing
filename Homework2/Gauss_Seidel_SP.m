function [x, errors] = Gauss_Seidel_SP(A,b, k)
% Implementation of Gauss-Seidel method for solving a system of equations of the
% form Ax = b. Will run for k iterations and return the solution as
% calculated and a list with error for each iteration.

    % Split the matrix in to lower and upper triangles
    L = tril(A);
    U = A-L;
    
    % allocate memory
    x = zeros(length(b),1);
    errors = zeros(1,k);
    
    for ii = 1:k
        x = inv(L) * (b - U*x); % calculate next x

        errors(ii) = norm(A*x - b); % find error
    end



end