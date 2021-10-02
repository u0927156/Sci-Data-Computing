function [x, errors] = Jacobi_SP(A,b, k)
% Implementation of Jacobi method for solving a system of equations of the
% form Ax = b. Will run for k iterations and return the solution as
% calculated and a list with error for each iteration.

    % Decompose the matrix into the two portions we need
    D = diag(diag(A));
    LU = A-D;

    % allocate the results
    x = zeros(length(b),1); % x_0 is all zeros
    errors = zeros(1,k);
    
    
    for ii = 1:k
        x = inv(D) * (b - LU *x); % calculate the next x

        errors(ii) = norm(A*x - b); % calculate error
    end



end