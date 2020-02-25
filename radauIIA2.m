function yi = radauIIA2(ti,yi,h,f,df)
tol=1e-14; maxit=2*100;
half=length(yi);
Phi = @(k) [k(1:half) - f(ti + h/3, yi + h*(5*k(1:half)/12 - k(half+1:end)/12));
            k(half+1:end) - f(ti + h, yi + h*(3*k(1:half)/4 + k(half+1:end)/4)) ];
dPhi = @(k) [ 1 - (5*h/12)*df(ti + h/3,yi + h*(5*k(1:half)/12 - k(half+1:end)/12)), h/12*df(ti + h/3,yi + h*(5*k(1:half)/12 - k(half+1:end)/12)); -(3*h/4)*df(ti + h,yi + h*(3*k(1:half)/4 + k(half+1:end)/4)), 1- (h/4)*df(ti + h, yi + h*(3*k(1:half)/4 + k(half+1:end)/4)) ];
[k, it] = newton(Phi,dPhi,[f(ti,yi);f(ti,yi)], tol, maxit);
yi = yi + h*(k(1:half)*3/4 + k(half+1:end)*1/4);
end %function
