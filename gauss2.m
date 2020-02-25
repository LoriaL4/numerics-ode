function yi = gauss2(ti,yi,h,f,df)
tol=1e-14; maxit=2*100;
half=length(yi);
Phi = @(k) [k(1:half) - f(ti + ((1/2) - (sqrt(3)/6))*h, yi + h*(k(1:half)/4 + ((1/4) - sqrt(3)/6)*k(half+1:end)));
            k(half+1:end) - f(ti + ((1/2) + (sqrt(3)/6))*h, yi + h*(((1/4) + sqrt(3)/6)*k(1:half) + k(half+1:end)/4)) ];
dPhi = @(k) [1-(h/4)*df(ti+((1/2)-(sqrt(3)/6))*h,yi+h*(k(1:half)/4 + ((1/4)-sqrt(3)/6)*k(half+1:end))) ,-h*((1/4)-sqrt(3)/6)*df(ti+((1/2)-(sqrt(3)/6))*h,yi+h*(k(1:half)/4 + ((1/4)-sqrt(3)/6)*k(half+1:end)));
             -h*((1/4)+sqrt(3)/6)*df(ti+((1/2)+(sqrt(3)/6))*h,yi+h*(((1/4)+sqrt(3)/6)*k(1:half) +k(half+1:end)/4)),1-(h/4)*df(ti+((1/2)+(sqrt(3)/6))*h,yi+h*(((1/4)+sqrt(3)/6)*k(1:half) + k(half+1:end)/4))
            ];
[k, it] = newton(Phi,dPhi,[f(ti,yi);f(ti,yi)], tol, maxit);
yi = yi + h*(k(1:half)/2 + k(half+1:end)/2);
end %function
