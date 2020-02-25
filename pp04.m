% pp04.m: 
% Gloria Doci
% Loesung der Wärmeleitungsgleichung
%            u_t(x,t)-u_xx(x,t) + u(x,t)=0, x in (0,1), t > 0
%                               u_x(x,t)=0, x in {0,1}, t > 0
% mit finiten Differenzen und konstanter Schrittweite
% Analytische Lösung u(x,t)= cn*cos(n*pi*x)*e^(t+n^2*pi^2*t) 

% Problembeschreibung

f = @(x,t) 0;

uex = @(x,t) cos(pi*x)*e^(-pi^2*t-t); % c1=1 and cn=0 for n >=2

u0= @(x) cos(pi*x);

 
T=1;
for theta=[0,0.5,1]
% Schleife über Anzahl der Elemente
ee=[]; hh=[]; tt=[];
fprintf("THETA=%f\n",theta)
for n=2.^[1:5]

% Gitter
h=1/n; x=[0:h:1]';

A=sparse(n+1,n+1); 

% die Randbedingung
A(1,1)=1/h^2*(h^2+2);
A(1,2)=-2/h^2;

A(n+1,n)=-2/h^2;
A(n+1,n+1)=1/h^2*(h^2+2);

% innere Punkte
for i=2:n
  A(i,i-1)=-1/h^2;
  A(i,i)=(h^2+2)/h^2;
  A(i,i+1)=-1/h^2;
end % for i


E=speye(n+1);

uh=u0(x);
tau=h^2/2;
err=0;
uh=u0(x);

for t=tau:tau:T
  uh=(E+theta*tau*A)\(uh-tau*(1-theta)*A*uh);
 
% Fehlerberechnung
eh=norm(uh-uex(x,t),'inf');

err=max(err,eh);
% Ausgabe

end
hh = [hh,h];
tt = [tt,t];
ee=[ee,err];

fprintf('n= %d h=%1.3e \t |u-uh|=%1.3e',n,h,err)
pause;

end

p1=polyfit(log(hh),log(ee),1);
fprintf('rate O(h^p):  \t         %f \t   \n',p1(1));
end
