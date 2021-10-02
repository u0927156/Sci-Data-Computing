function outMatrix = WaterMatrix(DeltaX, a, n)
% Makes a matix that represents the system of equations that arises from
% the flow of water through two very different materials. 

    % Preallocate the matrix
    outMatrix = zeros(n,n);

    % find the midpoint in advance
    midpoint = ceil(n/2);
    for ii = 1:n
        
        % Put in values for each of the different sections of matrix
        if ii==1
            outMatrix(1, 1:2) = [-2 1];
        elseif ii < midpoint
            outMatrix(ii, (ii-1):(ii+1)) = [1 -2 1];
        elseif ii == midpoint
            outMatrix(ii, (ii-1):(ii+1)) = [1 -(1+a) a];
        elseif ii > midpoint && ii ~= n
            outMatrix(ii, (ii-1):(ii+1)) = [a -2*a a];
        else
            outMatrix(ii, (ii-1):ii) = [a -2*a];
        end
    end
    
    % scale matrix
    outMatrix = outMatrix./DeltaX^2;
end 