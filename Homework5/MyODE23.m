function [t, y, error] = MyODE23(f, t0, tEnd, y0, h)
    % Implementation of ODE23 Method for problem 5 of homework 5 CS 6210. 
    
    curr_y = y0;
    curr_t = t0;
    ys = [curr_y];
    times = [curr_t];
    errors = [];
    
    while curr_t < tEnd
        S1 = f(curr_t, curr_y);

        S2 = f(curr_t + h/2, curr_y + h/2 * S1);

        S3 = f(curr_t + 3*h/4, curr_y + 3 * h/ 4 * S2);

        next_t = curr_t + h;

        next_y = curr_y + h/9 * (2 * S1 + 3 * S2 + 4 * S3);


        S4 = f(next_t, next_y);
        curr_error = h/72 * (-5 * S1 + 6 * S2 + 8*S3 - 9 * S4);

        times = [times next_t];
        ys = [ys next_y];

        errors = [errors curr_error];
        curr_t = next_t;
        curr_y = next_y;

    end
    t = times;
    y = ys;
    error = errors;
end