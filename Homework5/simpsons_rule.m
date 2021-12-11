function Area = simpsons_rule(F, a, b, N)
% Implementation of simpsons rule. F is a function handle, a and b are the
% starting and end points. N is the number of points to evaluate on. 
    width = (b-a)/N;

    points = width:width:b;

    Area = sum(F(points) .* width);

end