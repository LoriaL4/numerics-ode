% pp02.m
% Gloria Doci
% Lösung zu Programmierpraktikum 2 zur VL Numerik Gewoehnlicher Differentialgleichungen im WiSe 2018/19
source('ab.m')
fprintf('start p2.m \n')
%%%%%% Modellprobleme %%%%%
ex=2; % Beispielauswahl
if ex==1
  %  % (a) y' = (1+y^2), y(0)=0, T=1
  f = @(t,y) (1+y.^2); y0=0; T=1; yex = @(t)  tan(t);
  df = @(t,y) (2*y);
else
  % (b) y' = A y, A=[1 1; 0 1], y0=(1,1), T=1
  f = @(t,y) [1,1;0,1]*y; y0=[1;1]; T=1; yex = @(t) [(1+t).*exp(t); exp(t)];
  df = @(t,y) [1,1;0,1];
end

function [] = test_with_esv(start_points_esv,y0,T,yex,f,df,a,b,msv_name)
EE1h=[]; hh=[];

% loop over time steps
for nT=5*2.^(1:5)

h=T/nT; % mesh % T=1

% Startphase
yi = y0; ti = 0; tt = 0; yy = y0;
Eh = 0;

k=length(a)-1;

if start_points_esv==1
   func=@radauIIA1;
   esv_name='radauIIA1';
elseif start_points_esv==2
   func=@gauss1;
   esv_name='gauss1';
elseif start_points_esv==3
   func=@radauIIA2;
   esv_name='radauIIA2';
elseif start_points_esv==4
   func = @gauss2;
   esv_name='gauss2';
else error('invalid start points esv!')
end % if
for t = h:h:(k-1)*h
yi = feval(func,ti,yi,h,f,df);
ti = ti + h;
yy = [yy,yi]; 
tt = [tt,ti];
end

for j = 1:columns(yy)
Eh = max(Eh, norm(yex(tt(j))-yy(:,j)));
end %for

% jetzt weiter mit MSV
for t = k*h:h:T
[yy,tt] = msv_step(tt,yy,h,f,df,a,b);

% only one value added at the end of yy for each iteration
Eh = max(Eh, norm(yex(tt(end))-yy(:,end)));
end
fprintf(' h=%f: \t Eh_%s mit %s = %e \n',h,msv_name,esv_name,Eh)
EE1h=[EE1h,Eh]; hh=[hh,h];

end % for
p1=polyfit(log(hh),log(EE1h),1);
fprintf('rate O(h^p):  \t         %f \t \n',p1(1))
endfunction


fprintf('starting tests with Adams-Moulton k=2,p=3\n')
for i=1:4
test_with_esv(i,y0,T,yex,f,df,a_moulton_k2,b_moulton_k2,'moulton_k2');
end
fprintf('Press any key to continue with other tests\n')
pause;

fprintf('starting tests with Adams-Moulton k=3,p=4\n')
for i=1:4
test_with_esv(i,y0,T,yex,f,df,a_moulton_k3,b_moulton_k3,'moulton_k3');
end
fprintf('Press any key to continue with other tests\n')
pause;

fprintf('starting tests with Adams-Bashforth k=3,p=3\n')
for i=1:4
test_with_esv(i,y0,T,yex,f,df,a_bashforth_k3,b_bashforth_k3,'bashforth_k3');
end
fprintf('Press any key to continue with other tests\n')
pause;

fprintf('starting tests with BDF k=2,p=2\n')
for i=1:4
test_with_esv(i,y0,T,yex,f,df,a_bdf_k2,b_bdf_k2,'bdf_k2');
end

fprintf('Done\n')
