% pp03.m: 
% Gloria Doci
% Loesung des Poisson Problems
%            -u_xx(x,y)-u_yy(x,y) = f(x,y),  (x,y) in (0,1)^2
%               u(x)  = 0,                  (x,y) im Rand
% mit finiten Differenzen und konstanter Schrittweite

% Problembeschreibung

f = @(x,y) (2*pi^2)*sin(pi*x)*sin(pi*y); % f = @(x,y) 2*pi^2*sin(pi*x)*sin(pi*y);

uex = @(x,y) sin(pi*x)*sin(pi*y);

% Insgesamt (n+1)^2 Pukte
% Zeilenweise Numerierung, die von 1 anfängt
point = @(x,y,n) n*x+1+(n+1)*(y/(1/n)); 

findi = @(p,n) mod(p-1,n+1)*(1/n);

findj= @(p,n) idivide(p-1,n+1,"fix")*(1/n);  

% Schleife über Anzahl der Elemente
ee=[]; hh=[];
for n=2.^[1:5]

% Gitter
h=1/n; x=[0:h:1]'; y=[0:h:1]';

% Randbedingung bei (x,0),(x,1),(0,y),(1,y), d.h. alle außer inneren Punkten haben b(i)=0
b=zeros((n+1)^2,1);

tic
A=speye((n+1)^2); % Randbedingung bei (x,0),(x,1),(0,y),(1,y), d.h. alle außer inneren Punkten haben A(i,i)=1
 
% innere Punkte
for i=x(2:end-1)'
    for j=y(2:end-1)'

% Innere Punkte
  p=point(i,j,n);
  A(p,p-(n+1)) = -1/h^2; 
  A(p,p+(n+1)) = -1/h^2; 
  A(p,p-1)     = -1/h^2;
  A(p,p+1)     = -1/h^2;
  A(p,p)       = 4/h^2;
b(p) = f(i,j);
end % for j
end % for i

%{ 
Randpunkten
% Randpunkte (x,0), (x,1)
for i=x'
  pdown=point(i,0,n);
  pup=point(i,1,n);
  A(pdown,pdown)=1;
  A(pup,pup)=1;
end % for x

% Randpunkte (y,0), (y,1)
for j=y'
  pleft=point(0,j,n);
  pright=point(1,j,n);
  A(pleft,pleft) = 1;
  A(pright,pright) = 1;
end % for y
%}
elapsed_FD=toc;

tic
% Solve

tol=1e-14; maxit=2*100; 
[uh, ~] = pcg(A,b,tol,maxit,@(x) x);

elapsed_Solve=toc;

e=0;
% Fehlerberechnung
for i=x'
    for j=y'
    e=max(e,norm(uh(point(i,j,n))-uex(i,j),'inf'));
end
end

hh = [hh,h];
ee=[ee,e];
% Ausgabe
fprintf('n= %d h=%1.3e \t |u-uh|=%1.3e \t elapsed_FD_aufstellen=%f sec \t elapsed_LGS_lösen=%f sec \n',n,h,e,elapsed_FD,elapsed_Solve)


% plot
mesh(x,y,reshape(uh,n+1,n+1))
xlabel ("x");
ylabel ("y");
zlabel ("uh");
title (sprintf("uh for n=%d \t h=%f \t e=%f",n,h,e));
pause;
end
p1=polyfit(log(hh),log(ee),1);
fprintf('rate O(h^p):  \t         %f \t   \n',p1(1));
%[xx, yy] = meshgrid (x, y);
%mesh(x,y,uex(xx,yy))
%pause;
% Konvergenzplot
loglog(hh,0.02*hh.^2,'b--',hh,ee,'r*');
legend('0.02 h^2','|u-uh|','location','southeast')
