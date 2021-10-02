function [x, errors] = SOR_SP(A,b, omega, k)
% Implementation of successive over relaxation method for solving a system 
% of equations of the form Ax = b. Will run for k iterations and return 
% the solution as calculated and a list with error for each iteration.


    % split the matrix in to diagonal and lower triangle
    D = diag(diag(A));
    E = tril(A);
    
    % allocate memory
    x = zeros(length(b),1);
    errors = zeros(1,k);
    
    r = b - A*x;  % calculate initial residual
    for ii = 1:k
        x = x + omega * ((1-omega)*D + omega * E) \ r; % calculate next x

        r = b - A*x; % calculate next residual
        errors(ii) = norm(A*x - b); % calculate error
    end



end