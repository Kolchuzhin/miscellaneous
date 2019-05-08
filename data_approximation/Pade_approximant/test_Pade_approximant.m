%=========================================================================%
% PURPOSE:
%
%          scripts:  test_Pade_approximant.m
%                         Pade_approximant.m
%-------------------------------------------------------------------------%
% Kolchuzhin V.A., Aschheim, 01.03.2016
% <vladimir.kolchuzhin@ieee.org>
%=========================================================================%
% Some behavioral models cannot be presented efficiently by a Taylor series

% 1. singularity on the real axis: 1/x
% 2. singularities in the complex plane: 1/(1+x^2)

%-------------------------------------------------------------------------%
clear;
%-------------------------------------------------------------------------%
x0=1.0; y0=1.0/x0; % expansion point

X=[0.1:0.05:1.95]';
Y=1./(X);
%-------------------------------------------------------------------------%
% Taylor series expansion for the function 1/x about the point x0=1.0:
Tk=[1 -1 1 -1 1 -1]; 

k=numel(Tk)-1;
T=Tk(1); for i=1:k T=T+Tk(i+1).*(X-x0).^i; end
%-------------------------------------------------------------------------%
% Pade approximant
m=0; 
n=1;

[Pm,Qn] = Pade_approximant(Tk,m,n);

P=Pm(1); Q=Qn(1);
for i=1:m P=P+Pm(i+1).*(X-x0).^i; end
for i=1:n Q=Q+Qn(i+1).*(X-x0).^i; end
R=P./Q;
%-------------------------------------------------------------------------%
%
hold on;
grid on;
xlabel('x');
ylabel('f(x)');

plot(X,T,'-b','LineWidth',2);
plot(X,R,'-r','LineWidth',2);
plot(X,Y,'o','MarkerFaceColor','c');
plot(x0,y0,'or','MarkerFaceColor','y');

legend(['Taylor series: ',int2str(k),'order'],['Pade approximant: [',int2str(m),'/',int2str(n),']'],'data sampling: 1/x','expansion point');
%=========================================================================%
