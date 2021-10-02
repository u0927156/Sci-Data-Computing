function v_out = Vcycle_2d(v_in,f_in,nu_down,nu_up)
%
% assume M=N in the code, i.e. an equi-spaced square domain
%
[M N] = size(v_in);  M = M-1;  N = N-1;
if (M ~= N) 
   disp('M neq N!')
   v_out =  0*v_in;
   return
end

num_levels = log2(M);      % determine number of levels
ib = zeros(num_levels,1);  % begin index for each level 
ie = zeros(num_levels,1);  % end index for each level

ib(1) = 1; ie(1) = (M+1)^2;
len_arr = (M+1)^2;
cur_M = M;

for j = 2:num_levels
    cur_M = cur_M/2;
    ib(j) = ie(j-1)+1;
    ie(j) = ie(j-1)+(cur_M+1)^2;
    len_arr = len_arr+(cur_M+1)^2;
end

% initialize grid

v = zeros(len_arr,1); 
f = zeros(len_arr,1);
v(ib(1):ie(1)) = v_in(:);
f(ib(1):ie(1)) = f_in(:);

% V-cycle code

% traverse down grids

for j = 1:num_levels-1   

    ind_f = ib(j):ie(j);     % fine grid index range 
    ind_c = ib(j+1):ie(j+1); % coarse grid index range

    v(ind_f) = relax_2d(v(ind_f),f(ind_f),nu_down);
    f(ind_c) = restrict_2d(f(ind_f)-apply_operator_2d(v(ind_f)));

end  

% solve exactly on coarse grid 

h = 1/2; hsq = h*h; 
v(ib(num_levels):ie(num_levels)) = 0;
v(ib(num_levels)+4) = (hsq/4)*f(ib(num_levels)+4);

% traverse up grids
disp('going up')

for j = num_levels-1:-1:1   

    ind_f = ib(j):ie(j);     % fine grid index range 
    ind_c = ib(j+1):ie(j+1); % coarse grid index range

    v(ind_f) = v(ind_f) + interpolate_2d(v(ind_c));
    v(ind_c) = 0;
    v(ind_f) = relax_2d(v(ind_f),f(ind_f),nu_up);
    
end    

v_out = reshape(v(ib(1):ie(1)),M+1,M+1);
