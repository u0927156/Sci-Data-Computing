function y = interpolate_2d(v)
%
%
%
M = sqrt(length(v))-1;  
vm = reshape(v,M+1,M+1);
Mym = 2*M; ym = zeros(Mym+1);

% ym center points 

ym(3:2:Mym-1,3:2:Mym-1) = vm(2:M,2:M);

% ym SE, NE, SW & NW points

ym(2:2:Mym,2:2:Mym) = ( vm(1:M,  1:M) + vm(2:M+1,  1:M) + ...
                        vm(1:M,2:M+1) + vm(2:M+1,2:M+1) )/4;
% ym N & S points                    

ym(3:2:Mym-1,2:2:Mym) = ( vm(2:M,1:M) + vm(2:M,2:M+1) )/2;

% ym E & W points
                    
ym(2:2:Mym,3:2:Mym-1) = ( vm(1:M,2:M) + vm(2:M+1,2:M) )/2;

% put ym in y 

y = ym(:);



