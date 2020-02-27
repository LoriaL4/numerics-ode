function [k, it] = newton(G,dG,k0,tol,maxit)
% fprintf("start newton\n")
    k=k0;
    i=1;
    while i < maxit
        
        [delta, ~] = pcg(dG(k),-G(k),tol,maxit,@(x) x);
        % delta = -dG(k)\G(k);
        kii = k + delta;
        i = i + 1;  
        
        if norm(kii - k) <= tol 
           k = kii;
                      
           break;                 
        else        
         k = kii; 
        end % end if                     
    end % end while
    it = i;        
end % function
