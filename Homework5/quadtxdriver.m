fprintf('  tol            Q                 fcount      err       ratio \n')

for  k = 1:12
    %kk = 6
    tol = 10^(-k);
 %humps   Qexact=29.85832539549867;
    [Q,fcount] = quadtx(@burger,0,1,tol);
 %   err=Q-Qexact;
 %   ratio = err/tol;
 fprintf('%8.0e %21.14f %7d  \n',tol,Q,fcount)
 %   fprintf('%8.0e %21.14f %7d %13.3e %9.3f \n',tol,Q,fcount,err,ratio)
end