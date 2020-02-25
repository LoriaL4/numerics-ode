% pp0.m
% Loesungsvorschlag zu Programmierpraktikum 0 zur
% VL Numerik Gewoehnlicher Differentialgleichungen
% im WiSe 2018/19

%%%%%% Modellprobleme %%%%%
ex=2; % Beispielauswahl
if ex==1
  %  % (a) y' = (1+y^2), y(0)=0, T=1
  f = @(t,y) (1+y.^2); y0=0; T=1; yex = @(t)  tan(t);
else
  % (b) y' = A y, A=[1 1; 0 1], y0=(1,1), T=1
  f = @(t,y) [1,1;0,1]*y; y0=[1;1]; T=1; yex = @(t) [(1+t).*exp(t); exp(t)];
end

%%%%%% Numerical Integration %%%%%

% store errors
EE1h=[]; EE1hw=[]; 
EE2h=[]; EE2hw=[]; 
EE3h=[]; EE3hw=[]; 
EE4h=[]; hh=[];

% loop over time steps
for nT=5*2.^(1:5)

h=T/nT; % mesh % T=1

% Initialize
ti=0; 
y1i=y0; % exp Euler 
y2i=y0; % verb Euler
y3i=y0; % Heun-3 
y4i=y0; % Runge-Kutta-4

% time stepping
E1h=0; E1hw=0; E2h=0; E2hw=0; E3h=0; E3hw=0; E4h=0;
for i=2:nT+1

    % Exp Euler 
    y1i = y1i + h*f(ti,y1i);

    % verb Euler
    k1  = f(ti,y2i);
    k2  = f(ti+h/2,y2i+h/2*k1); 
    y2i = y2i + h*k2;       

    % Heun-3
    k1  = f(ti,y3i);
    k2  = f(ti+h/3,y3i+h/3*k1); 
    k3  = f(ti+2*h/3,y3i+2*h/3*k2);
    y3i = y3i + h*(0.25*k1+0.75*k3);

    % Runge-Kutta-4
    k1  = f(ti,y4i);
    k2  = f(ti+h/2,y4i+h/2*k1); 
    k3  = f(ti+h/2,y4i+h/2*k2);
    k4  = f(ti+h,y4i+h*k3);
    y4i = y4i + h*(k1+2*k2+2*k3+k4)/6;

    % time-step
    ti  = ti + h;                     

    % Fehlerberechnung
    E1h  = max(E1h,norm(yex(ti)-y1i));
    E1hw = max(E1hw,norm(y2i-y1i));
    %
    E2h  = max(E2h,norm(yex(ti)-y2i));
    E2hw = max(E2hw,norm(y3i-y2i));
    %
    E3h  = max(E3h,norm(yex(ti)-y3i));
    E3hw = max(E3hw,norm(y4i-y3i));
    %
    E4h  = max(E4h,norm(yex(ti)-y4i));
   
end 

% Fehlerausgabe
fprintf(' h=%f: \t E1h=%e \t E1hw=%e \t E2h=%e \t E2hw=%e \t E3h=%e \t E3hw=%e \t E4h=%e\n',h,E1h,E1hw,E2h,E2hw,E3h,E3hw,E4h)

EE1h=[EE1h,E1h]; EE1hw=[EE1hw,E1hw]; 
EE2h=[EE2h,E2h]; EE2hw=[EE2hw,E2hw]; 
EE3h=[EE3h,E3h]; EE3hw=[EE3hw,E3hw]; 
EE4h=[EE4h,E4h]; hh=[hh,h];

end % loop over time steps
p1=polyfit(log(hh),log(EE1h),1); p1w=polyfit(log(hh),log(EE1hw),1);
p2=polyfit(log(hh),log(EE2h),1); p2w=polyfit(log(hh),log(EE2hw),1);
p3=polyfit(log(hh),log(EE3h),1); p3w=polyfit(log(hh),log(EE3hw),1);
p4=polyfit(log(hh),log(EE4h),1); 
fprintf('rate O(h^p):  \t         %f \t          %f \t         %f \t          %f \t         %f \t          %f \t         %f\n',p1(1),p1w(1),p2(1),p2w(1),p3(1),p3w(1),p4(1))

% Konvergenzplot
loglog(hh,hh,'b--',hh,EE1h,'b*',hh,EE1hw,'bo',...
       hh,hh.^2,'r--',hh,EE2h,'r*',hh,EE2hw,'ro',...
       hh,hh.^3,'c--',hh,EE3h,'c.',hh,EE3hw,'co',...
       hh,hh.^4,'k--',hh,EE4h,'k.');
title('Konvergenzverhalten')
legend('h','E1h','E1hw','h^2','E1h','E1hw','h^3','E1h','E1hw','h^4','E1h','location','northwest')

% Diskussion:
% Man sieht, dass die Raten genau den theoretischen Raten entsprechen.
% Die geschaetzten Fehler Ehw stimmen ebenfalls sehr gut mit den tatsaechlichen
% Fehler ueberein. 

