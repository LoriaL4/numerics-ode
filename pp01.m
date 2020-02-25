% pp01.m
% Gloria Doci
% Lösung zu Programmierpraktikum 1 zur VL Numerik Gewoehnlicher Differentialgleichungen im WiSe 2018/19
fprintf('start\n')
%%%%%% Modellprobleme %%%%%
ex=1; % Beispielauswahl
fprintf('you chose ex=%d\n',ex)
if ex==1
  %  % (a) y' = (1+y^2), y(0)=0, T=1
  f = @(t,y) (1+y.^2); y0=0; T=1; yex = @(t)  tan(t);
  df = @(t,y) (2*y);
else
  % (b) y' = A y, A=[1 1; 0 1], y0=(1,1), T=1
  f = @(t,y) [1,1;0,1]*y; y0=[1;1]; T=1; yex = @(t) [(1+t).*exp(t); exp(t)];
  df = @(t,y) [1,1;0,1];
end
 

%{
f = @(y) y^3-2*y-5;
df = @(y) 3*y^2-2;
[k, it] = newton(f, df, 0.9, 0.001, 100)
F=@(Y)   [2*Y(1);
          Y(2)^3-2*Y(2)-5];
dF=@(Y) [2 0; 0 3*Y(2)^2-2];
[k, it] = newton(F, dF, [0; 0], 0.001, 100)
%}

% P3 a)
h = 2^-1;
G1 = @(k) k-sin(1 + h*k);
dG1 = @(k) 1-cos(1 + h*k)*h;


%[k, it] = newton(G1, dG1, G1(0), 0.001, 100)

% P3 b)
G2 = @(k)  [k(1) - sin(h*k(1)) + cos(h*k(2));
            k(2) - cos(h*k(1)) - sin(h*k(2))];
dG2 = @(k) [1 - cos(h*k(1))*h, -sin(h*k(2))*h ; sin(h*k(1))*h, 1-cos(h*k(2))*h];

% [k, it] = newton(G2, dG2, [0;0], 0.001, 1000)

% P3 c)


G3 = @(k) k-A*(yi + h*k);
%G3([1; 0; 0]);
% toDo dG3

% store errors
EE1h=[]; EE2h=[]; EE3h=[]; EE4h=[]; 
hh=[];
fprintf('\t  \t radauIIA1 \t \t gauss1 \t \t radauIIA2 \t \t gauss2\n')
% loop over time steps
for nT=5*2.^(1:5)

h=T/nT; % mesh % T=1


x_vector_1=[]; x_vector_2=[]; x_vector_3=[]; x_vector_4=[];
y_vector_1=[]; y_vector_2=[]; y_vector_3=[]; y_vector_4=[];
y_ex=[];

yi_1=y0; yi_2=y0; yi_3=y0; yi_4=y0; ti=0;
% time stepping

E1h=0; E2h=0; E3h=0; E4h=0; 
for i=2:nT+1
    yi_1=radauIIA1(ti,yi_1,h,f,df);
    yi_2=gauss1(ti,yi_2,h,f,df);
    yi_3=radauIIA2(ti,yi_3,h,f,df);
    yi_4=gauss2(ti,yi_4,h,f,df);

    x_vector_1 = [x_vector_1, ti];
    y_vector_1 = [y_vector_1, yi_1];
    x_vector_2 = [x_vector_2, ti];
    y_vector_2 = [y_vector_2, yi_2];
    x_vector_3 = [x_vector_3, ti];
    y_vector_3 = [y_vector_3, yi_3];
    x_vector_4 = [x_vector_4, ti];
    y_vector_4 = [y_vector_4, yi_4];
    % y_ex = [y_ex, yex(ti)];
    ti = ti + h;
    % yi
    E1h  = max(E1h,norm(yex(ti)-yi_1));
    E2h  = max(E2h,norm(yex(ti)-yi_2));
    E3h  = max(E3h,norm(yex(ti)-yi_3));
    E4h  = max(E4h,norm(yex(ti)-yi_4));
end %for
% Fehlerausgabe

fprintf(' h=%f: \t E1h=%e \t E2h=%e \t E3h=%e \t E4h=%e \n',h,E1h,E2h,E3h,E4h)
EE1h=[EE1h,E1h]; EE2h=[EE2h,E2h]; EE3h=[EE3h,E3h]; EE4h=[EE4h,E4h]; hh=[hh,h];

end % for time steps
p1=polyfit(log(hh),log(EE1h),1);
p2=polyfit(log(hh),log(EE2h),1);
p3=polyfit(log(hh),log(EE3h),1);
p4=polyfit(log(hh),log(EE4h),1);
fprintf('rate O(h^p):  \t         %f \t         %f \t         %f \t         %f\n',p1(1),p2(1),p3(1),p4(1))
% plot(x_vector, y_vector, 'r-')
% hold on;
%plot(x_vector, y_ex, 'b-')
