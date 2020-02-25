function yi = radauIIA1(ti,yi,h,f,df)
tol=1e-14; maxit=2*100;
Phi = @(k) k - f(ti + h, yi + h*k);
dPhi = @(k) diag(ones(length(yi),1)) - h*df(ti + h, yi + h*k); % diag(ones(length(k0),1))
[k, it] = newton(Phi,dPhi,f(ti,yi), tol, maxit);
yi = yi + h*k;
end %function
