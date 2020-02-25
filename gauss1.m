function yi = gauss1(ti,yi,h,f,df)
tol=1e-14; maxit=2*100;
Phi = @(k) k - f(ti + h/2, yi + h*k/2);
dPhi = @(k) eye(length(yi)) - (h/2)*df(ti + h/2, yi + h*k/2);
[k, it] = newton(Phi, dPhi, f(ti,yi), tol, maxit);
yi = yi + h*k;
end %function
