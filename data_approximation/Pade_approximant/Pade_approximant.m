function [P,Q] = Pade_approximant(T,m,n)
%=========================================================================%
% PURPOSE:
%           univariate Pade approximant: Tk -> Pm/Qn
%
% INPUT:
%           T - taylor polynomial: [t0 t1 t2 ...]
%           m - max degree of numerator P
%           n - max degree of denominator Q
%
%           (m+n)<=k - max degree of polynomial T
% OUTPUT:
%           P - array of coefficients: [p0 p1 p2 ...]
%           Q - array of coefficients: [q0 q1 q2 ...]
%-------------------------------------------------------------------------%
% written by Kolchuzhin V.A., LMGT, TU Chemnitz, 14:08 11.07.2003
% <vladimir.kolchuzhin@etit.tu-chemnitz.de>
%=========================================================================%
% ver. 1.0 tested
%-------------------------------------------------------------------------%
if nargin==0 % self-test data
%------------
T=[1 2 3 4 5]; % k=4
m=2; n=2;

[P,Q] = Pade_approximant(T,m,n);
%------------
% P=[1  0 0];
% Q=[1 -2 1];
%------------
end
%-------------------------------------------------------------------------%
% m+n terms of polynomial T
k=m+n;
for i=1:k+1 A(i)=T(i); end
%-----------------------
% assembling K
K=zeros(m+n+1,m+n+1);
%==== 1 step ====
for i=1:m+1 K(i,i)=1; end
%==== 2 step ====
di=1;
for j=m+2:m+n+1
    for i=1+di:m+n+1
        K(i,j)=-A(i-di);
    end 
    di=di+1;
end
%-----------------------
% solution
pq=inv(K)*A';
%-----------------------
% substitution P and Q 
P=zeros(m+1,1);
Q=zeros(n+1,1); Q(1)=1;
for i=1:m+1 P(i)=pq(i);   end
for i=2:n+1 Q(i)=pq(i+m); end
%-------------------------------------------------------------------------%
