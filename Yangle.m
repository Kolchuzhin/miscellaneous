% A Matlab script to calculate silicon Young'g modulus and Poisson's ration in
% any crystal direction. For explation please see the companion tutorial
% at http://www.kaajakari.net/~ville/research/tutorials/elasticity_tutorial.pdf
%
% You may modify this script to your own use and you may distribute this script 
% unmodified. 
% 
% Author: Ville Kaajakari
% Date: 2.4.2003

c1111=1.66;
c1122=0.64;
c2323=0.80;

%Non-zero componenents of Si stiffness tensor in [100] coordinates
Cc(1,1,1,1)=c1111;
Cc(2,2,2,2)=c1111;
Cc(3,3,3,3)=c1111;

Cc(1,1,2,2)=c1122;
Cc(1,1,3,3)=c1122;
Cc(2,2,1,1)=c1122;
Cc(2,2,3,3)=c1122;
Cc(3,3,1,1)=c1122;
Cc(3,3,2,2)=c1122;

Cc(1,3,1,3)=c2323;
Cc(3,1,1,3)=c2323;
Cc(1,3,3,1)=c2323;
Cc(3,1,3,1)=c2323;

Cc(2,3,2,3)=c2323;
Cc(3,2,2,3)=c2323;
Cc(2,3,3,2)=c2323;
Cc(3,2,3,2)=c2323;

Cc(1,2,1,2)=c2323;
Cc(2,1,1,2)=c2323;
Cc(1,2,2,1)=c2323;
Cc(2,1,2,1)=c2323;

N=101;
theta=linspace(0,pi/2,N);
v_angle1=zeros(1,N);
v_angle2=zeros(1,N);

for kk=1:N,
    
    % Rotation matrix
    Q=[
        [cos(theta(kk)) sin(theta(kk)) 0];
        [-sin(theta(kk)) cos(theta(kk)) 0];
        [0 0 1];
    ];
    
    %Carculate stiffness tensor in rotated coordinates
    Crot=zeros(3,3,3,3);
    for i=1:3,
        for j=1:3,
            for k=1:3,
                for l=1:3,
                    for p=1:3,
                        for q=1:3,
                            for r=1:3,
                                for s=1:3,
                                    Crot(i,j,k,l)=Crot(i,j,k,l)+Q(p,i)*Q(q,j)*Q(r,k)*Q(s,l)*Cc(p,q,r,s);     
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    Cr=[
        [Crot(1,1,1,1) Crot(1,1,2,2) Crot(1,1,3,3) Crot(1,1,2,3) Crot(1,1,1,3) Crot(1,1,1,2)];
        [Crot(2,2,1,1) Crot(2,2,2,2) Crot(2,2,3,3) Crot(2,2,2,3) Crot(2,2,1,3) Crot(2,2,1,2)];
        [Crot(3,3,1,1) Crot(3,3,2,2) Crot(3,3,3,3) Crot(3,3,2,3) Crot(3,3,1,3) Crot(3,3,1,2)];
        [Crot(3,2,1,1) Crot(3,2,2,2) Crot(3,2,3,3) Crot(2,3,2,3) Crot(2,3,1,2) Crot(2,3,1,2)];
        [Crot(3,1,1,1) Crot(3,1,2,2) Crot(3,1,3,3) Crot(3,1,2,3) Crot(1,3,1,3) Crot(3,1,1,2)];
        [Crot(2,1,1,1) Crot(2,1,2,2) Crot(2,1,3,3) Crot(2,1,2,3) Crot(2,1,1,3) Crot(1,2,1,2)];
    ];
    
    % By definitions
    S=inv(Cr);
    Y_angle(kk)=1/S(1,1);
    v_angle1(kk)=-S(1,2)/S(1,1); % Between two (100)-plane vectors
    v_angle2(kk)=-S(1,3)/S(1,1); % Between (100)-plane vector and [001]-vector
    
end

figure(1);
h=polar(theta,Y_angle);
ha=axis;
axis([0 ha(2) 0 ha(4)]);

figure(2);
polar(theta,v_angle2)
hold on;
polar(theta,v_angle1);
hold off;
ha=axis;
axis([0 ha(2) 0 ha(4)]);
