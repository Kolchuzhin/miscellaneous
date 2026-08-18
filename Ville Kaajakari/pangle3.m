% A Matlab script to calculate silicon piezoresistance coefficients in
% any crystal direction. This script was originally downloaded from
% http://www.kaajakari.net/PracticalMEMS
%
% You may modify this script to your own use and you may distribute this script 
% unmodified. 
% 
% Author: Ville Kaajakari
% Date: 2.4.2009

function pangle3

%100
pr=protated(0,0,'n');
% By definitions
p_l_100n=pr(1,1)
p_t_100n=pr(1,2)
pr=protated(0,0,'p');
% By definitions
p_l_100p=pr(1,1)
p_t_100p=pr(1,2)

%110
pr=protated(pi/4,0,'n');
% By definitions
p_l_110n=pr(1,1)
p_t_110n=pr(1,2)
pr=protated(pi/4,0,'p');
% By definitions
p_l_110p=pr(1,1)
p_t_110p=pr(1,2)

%111
pr=protated(-pi/4,-asin(1/sqrt(3)),'n');
% By definitions
p_l_111n=pr(1,1)
p_t_111n=pr(1,2)
pr=protated(-pi/4,-asin(1/sqrt(3)),'p');
% By definitions
p_l_111p=pr(1,1)
p_t_111p=pr(1,2)

%110
N=51;
theta=linspace(0,pi,N);
Y_angle=zeros(1,N);
v_angle=zeros(1,N);


for kk=1:N,
   
    pr=protated(theta(kk),0,'n');

    % By definitions
    p_l(kk)=pr(1,1);
    p_t(kk)=pr(1,2);

    if p_l(kk)<0,
        theta_l(kk)=theta(kk)+pi;
    else
        theta_l(kk)=theta(kk);
    end

    if p_t(kk)<0,
        theta_t(kk)=theta(kk)+pi;
    else
        theta_t(kk)=theta(kk);
    end
end
%break

h=polar(theta_l,abs(p_l)); hold on
h=polar(theta_t,abs(p_t)); hold off

% mmpolar.m function can be used for nicer plots. Try Googling!
%h=mmpolar(theta_l,abs(p_l),'b',theta_t,abs(p_t),'r');
%mmpolar('RTickValue',[20 40 60 80 100 120])
%set(h,'linewidth',1.5);

function pr=protated(theta,gamma,type)

if type=='n'
    % n-type (ok)
    p1111=-102.2;
    p1122=53.4;
    p2323=-13.6/2;
else
    % p-type
    p1111=6.6;
    p1122=-1.1;
    p2323=138.1/2;
end
%Non-zero componenents of Si stiffness tensor in [100] coordinates
pp(1,1,1,1)=p1111;
pp(2,2,2,2)=p1111;
pp(3,3,3,3)=p1111;

pp(1,1,2,2)=p1122;
pp(1,1,3,3)=p1122;
pp(2,2,1,1)=p1122;
pp(2,2,3,3)=p1122;
pp(3,3,1,1)=p1122;
pp(3,3,2,2)=p1122;

pp(1,3,1,3)=p2323;
pp(3,1,1,3)=p2323;
pp(1,3,3,1)=p2323;
pp(3,1,3,1)=p2323;

pp(2,3,2,3)=p2323;
pp(3,2,2,3)=p2323;
pp(2,3,3,2)=p2323;
pp(3,2,3,2)=p2323;

pp(1,2,1,2)=p2323;
pp(2,1,1,2)=p2323;
pp(1,2,2,1)=p2323;
pp(2,1,2,1)=p2323;

p0=[
    [pp(1,1,1,1) pp(1,1,2,2) pp(1,1,3,3) pp(1,1,2,3) pp(1,1,1,3) pp(1,1,1,2)];
    [pp(2,2,1,1) pp(2,2,2,2) pp(2,2,3,3) pp(2,2,2,3) pp(2,2,1,3) pp(2,2,1,2)];
    [pp(3,3,1,1) pp(3,3,2,2) pp(3,3,3,3) pp(3,3,2,3) pp(3,3,1,3) pp(3,3,1,2)];
    [pp(3,2,1,1) pp(3,2,2,2) pp(3,2,3,3) pp(2,3,2,3) pp(2,3,1,2) pp(2,3,1,2)];
    [pp(3,1,1,1) pp(3,1,2,2) pp(3,1,3,3) pp(3,1,2,3) pp(1,3,1,3) pp(3,1,1,2)];
    [pp(2,1,1,1) pp(2,1,2,2) pp(2,1,3,3) pp(2,1,2,3) pp(2,1,1,3) pp(1,2,1,2)];
    ];


% Rotation matrix
Q1=[
    [cos(theta) sin(theta) 0];
    [-sin(theta) cos(theta) 0];
    [0 0 1];
    ];
l1=cos(gamma); m1=0; n1=sin(gamma);
l2=0;             m2=1;             n2=0;
l3=-sin(gamma); m3=0; n3=cos(gamma);

Q2=[
    [l1 m1 n1];
    [l2 m2 n2];
    [l3 m3 n3];
];
Q=Q1*Q2;
%Carculate stiffness tensor in rotated coordinates
prot=zeros(3,3,3,3);
for i=1:3,
    for j=1:3,
        for k=1:3,
            for l=1:3,
                for p=1:3,
                    for q=1:3,
                        for r=1:3,
                            for s=1:3,
                                prot(i,j,k,l)=prot(i,j,k,l)+Q(p,i)*Q(q,j)*Q(r,k)*Q(s,l)*pp(p,q,r,s);
                            end
                        end
                    end
                end
            end
        end
    end
end

pr=[
    [prot(1,1,1,1) prot(1,1,2,2) prot(1,1,3,3) prot(1,1,2,3) prot(1,1,1,3) prot(1,1,1,2)];
    [prot(2,2,1,1) prot(2,2,2,2) prot(2,2,3,3) prot(2,2,2,3) prot(2,2,1,3) prot(2,2,1,2)];
    [prot(3,3,1,1) prot(3,3,2,2) prot(3,3,3,3) prot(3,3,2,3) prot(3,3,1,3) prot(3,3,1,2)];
    [prot(3,2,1,1) prot(3,2,2,2) prot(3,2,3,3) prot(2,3,2,3) prot(2,3,1,2) prot(2,3,1,2)];
    [prot(3,1,1,1) prot(3,1,2,2) prot(3,1,3,3) prot(3,1,2,3) prot(1,3,1,3) prot(3,1,1,2)];
    [prot(2,1,1,1) prot(2,1,2,2) prot(2,1,3,3) prot(2,1,2,3) prot(2,1,1,3) prot(1,2,1,2)];
    ];



