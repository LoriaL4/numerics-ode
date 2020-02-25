function [yy,tt] = msv_step(tt,yy,h,f,df,a,b)
kk = length(a);
rhs = 0;
lhs = 0;
for m = 1:1:length(tt)
    rhs = rhs + b(m)*f(tt(m),yy(:,m));
    lhs = lhs + a(m)*yy(:,m);
end % for

tiplusk = tt(end) + h;
tt  = [tt, tiplusk];

if b(kk) == 0
y = (-lhs + h*rhs)/a(kk); % a(kk) ist immer != 0
else
G = @(yiplusk) a(kk)*yiplusk + lhs - h*rhs - h*b(kk)*f(tt(end),yiplusk);  
dG = @(yiplusk) a(kk)*eye(length(yy(:, end))) - h*b(kk)*df(tt(end), yiplusk);

[y,it] = newton(G,dG,yy(:,end),1e-14,200);
end % if
yy = [yy, y];
yy=yy(:,2:end);
tt=tt(2:end);
end % function
