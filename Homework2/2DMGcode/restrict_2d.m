function y = restrict_2d(v)
%
% 2D full-weighting restriction operator
%
M = sqrt(length(v))-1;  
vm = reshape(v,M+1,M+1);
Mym = M/2; ym = zeros(Mym+1);

ind = 3:2:M-1;  
ym(2:Mym,2:Mym) = ( 4*vm(ind,ind) + ...
                    2*(vm(ind-1,ind)+vm(ind+1,ind)+ ...
                       vm(ind,ind-1)+vm(ind,ind+1)) + ...
                      (vm(ind-1,ind-1)+vm(ind-1,ind+1)+ ...
                       vm(ind+1,ind-1)+vm(ind+1,ind+1)) )/16;
% put ym in y 

y = ym(:);


