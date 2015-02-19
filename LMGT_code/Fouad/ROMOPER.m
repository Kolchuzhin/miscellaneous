function [F,dF,ddF]=ROMOPER(q,PCOE,ORDER,FUN,TYPE)
%=========================================================================%
% PURPOSE:
%
%
% INPUT:
%            q,PCOE,ORDER,FUN,TYPE
% OUTPUT:
%            F,dF,ddF
% REFERENCE:
%            http://archiv.tu-chemnitz.de/pub/2005/0136
%            Dr.-Ing. F.Bennini, Dissertation, TU Chemnitz, 2005
%            A.2 MATLAB-Funktion zur numerischen Ermittlung der 
%            Funktionswerte und der partiellen Ableitungen
%-------------------------------------------------------------------------%
% Tested by Kolchuzhin V.A., LMGT, TU Chemnitz, 24.01.2011 13:57
% rev. 24.01.2011 11:42
%=========================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   q:     Vektor(1:3)   --> Modalen Verschiebungen                       %
%   PCOE:  Vektor        --> Polynomkoeffizienten                         %
%   ORDER: VEKTOR(1:3)   --> Polynom-Ordnung                              %
%   FUN:   Ganzzahl      --> 1 ... Vollbesetzter Polynom                  %
%                            2 ... Pascal Polynom                         %
%                            3 ... Reduzierter Polynom                    %
%   TYPE:  Ganzzahl      --> 1 ... Normal                                 %
%                            2 ... Inverse                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========================================================================%
%global REGFAC
REGFAC = [1 1 1];   % !!!!

  W=0;
 dW=zeros(1,3);
ddW=zeros(3,3);
%------------------------------------ Polynomgliedern der erste Variable
for iord=0:ORDER(1)
    x(iord+1)=q(1)^(iord)*REGFAC(1)^(iord);
    %Erste Ableitung
    if iord==1
        dx(iord+1)=REGFAC(1);
    elseif iord-1<0
        dx(iord+1)=0;
    else
        dx(iord+1)=iord*q(1)^(iord-1)*REGFAC(1)^iord;
    end
    %Zweite Ableitung
    if iord==2
        ddx(iord+1)=2*REGFAC(1)^2;
    elseif iord-1<0 || iord-2<0
        ddx(iord+1)=0;
    else
        ddx(iord+1)=iord*(iord-1)*q(1)^(iord-2)*REGFAC(1)^iord;
    end
end
%------------------------------------ Polynomgliedern der zweite Variable
for iord=0:ORDER(2)
    y(iord+1)=q(2)^(iord)*REGFAC(2)^(iord);
    %Erste Ableitung
    if iord==1
        dy(iord+1)=REGFAC(2);
    elseif iord-1<0
        dy(iord+1)=0;
    else
        dy(iord+1)=iord*q(2)^(iord-1)*REGFAC(2)^iord;
    end
    %Zweite Ableitung
    if iord==2
        ddy(iord+1)=2*REGFAC(2)^2;
    elseif iord-1<0 || iord-2<0
        ddy(iord+1)=0;
    else
        ddy(iord+1)=iord*(iord-1)*q(2)^(iord-2)*REGFAC(2)^iord;
    end
end
%------------------------------------ Polynomgliedern der dritte Variable
for iord=0:ORDER(3)
    z(iord+1)=q(3)^(iord)*REGFAC(3)^(iord);
    %Erste Ableitung
    if iord==1
        dz(iord+1)=REGFAC(3);
    elseif iord-1<0
        dz(iord+1)=0;
    else
        dz(iord+1)=iord*q(3)^(iord-1)*REGFAC(3)^iord;
    end
    %Zweite Ableitung
    if iord==2
        ddz(iord+1)=2*REGFAC(3)^2;
    elseif iord-1<0 || iord-2<0
        ddz(iord+1)=0;
    else
        ddz(iord+1)=iord*(iord-1)*q(3)^(iord-2)*REGFAC(3)^iord;
    end
end
%------------------------------------ Polynomaufbau mit drei Variablen
if FUN==1                           % vollbesetztes Polynom
    count=1;
    for kord=ORDER(3):-1:0          %z
        for jord=ORDER(2):-1:0      %y
            for iord=ORDER(1):-1:0  %x
                W=W+PCOE(count)*x(iord+1)*y(jord+1)*z(kord+1);
                dW(1)=dW(1)+PCOE(count)*dx(iord+1)*y(jord+1)*z(kord+1);
                dW(2)=dW(2)+PCOE(count)*x(iord+1)*dy(jord+1)*z(kord+1);
                dW(3)=dW(3)+PCOE(count)*x(iord+1)*y(jord+1)*dz(kord+1);
                ddW(1,1)=ddW(1,1)+PCOE(count)*ddx(iord+1)*y(jord+1)*z(kord+1);
                ddW(1,2)=ddW(1,2)+PCOE(count)*dx(iord+1)*dy(jord+1)*z(kord+1);
                ddW(1,3)=ddW(1,3)+PCOE(count)*dx(iord+1)*y(jord+1)*dz(kord+1);
                ddW(2,1)=ddW(2,1)+PCOE(count)*dx(iord+1)*dy(jord+1)*z(kord+1);
                ddW(2,2)=ddW(2,2)+PCOE(count)*x(iord+1)*ddy(jord+1)*z(kord+1);
                ddW(2,3)=ddW(2,3)+PCOE(count)*x(iord+1)*dy(jord+1)*dz(kord+1);
                ddW(3,1)=ddW(3,1)+PCOE(count)*dx(iord+1)*y(jord+1)*dz(kord+1);
                ddW(3,2)=ddW(3,2)+PCOE(count)*x(iord+1)*dy(jord+1)*dz(kord+1);
                ddW(3,3)=ddW(3,3)+PCOE(count)*x(iord+1)*y(jord+1)*ddz(kord+1);
                count=count+1;
            end
        end
    end
elseif FUN==2                       % Pascal Polynom
    count=1;
    Pord=ORDER;
    for kord=0:Pord(3)
        Pord(1)=ORDER(1)-kord;
        Pord(2)=ORDER(2);
        for jord=0:Pord(2)
            for iord=0:Pord(1)
                W=W+PCOE(count)*x(iord+1)*y(jord+1)*z(kord+1);
                dW(1)=dW(1)+PCOE(count)*dx(iord+1)*y(jord+1)*z(kord+1);
                dW(2)=dW(2)+PCOE(count)*x(iord+1)*dy(jord+1)*z(kord+1);
                dW(3)=dW(3)+PCOE(count)*x(iord+1)*y(jord+1)*dz(kord+1);
                ddW(1,1)=ddW(1,1)+PCOE(count)*ddx(iord+1)*y(jord+1)*z(kord+1);
                ddW(2,2)=ddW(2,2)+PCOE(count)*x(iord+1)*ddy(jord+1)*z(kord+1);
                ddW(3,3)=ddW(3,3)+PCOE(count)*x(iord+1)*y(jord+1)*ddz(kord+1);
                if iord==1 && jord>1
                    ddW(1,2)=ddW(1,2)+PCOE(count)*REGFAC(1)...
                    *(jord*(REGFAC(2)^jord)*(q(2)^(jord-1)))*((REGFAC(3)^kord)*(q(3)^kord));
                elseif jord==1 && iord>1
                    ddW(1,2)=ddW(1,2)+PCOE(count)*REGFAC(2)...
                    *(iord*(REGFAC(1)^iord)*(q(1)^(iord-1)))*((REGFAC(3)^kord)*(q(3)^kord));
                elseif jord==1 && iord==1
                    ddW(1,2)=ddW(1,2)+PCOE(count)*REGFAC(2)...
                    *REGFAC(1)*((REGFAC(3)^kord)*(q(3)^kord));
                elseif iord-1<0 || jord-1<0
                    ddW(1,2)=ddW(1,2)+0;
                else
                    ddW(1,2)=ddW(1,2)+PCOE(count)*iord*(REGFAC(1)^iord)...
                    *(q(1)^(iord-1))*(jord*(REGFAC(2)^jord)*(q(2)^(jord-1)))...
                *((REGFAC(3)^kord)*(q(3)^kord));
                end
                ddW(2,1)=ddW(1,2);
                if iord==1 && kord>1
                    ddW(1,3)=ddW(1,3)+PCOE(count)*REGFAC(1)*(kord*(REGFAC(3)^kord)...
                    *(q(3)^(kord-1)))*((REGFAC(2)^jord)*(q(2)^jord));
                elseif kord==1 && iord>1
                    ddW(1,3)=ddW(1,3)+PCOE(count)*REGFAC(3)*(iord*(REGFAC(1)^iord)...
                    *(q(1)^(iord-1)))*((REGFAC(2)^jord)*(q(2)^jord));
                elseif kord==1 && iord==1
                    ddW(1,3)=ddW(1,3)+PCOE(count)*REGFAC(3)*REGFAC(1)*((REGFAC(2)^jord)...
                    *(q(2)^jord));
                elseif iord-1<0 || kord-1<0
                    ddW(1,3)=ddW(1,3)+0;
                else
                    ddW(1,3)=ddW(1,3)+PCOE(count)*iord*(REGFAC(1)^iord)*(q(1)^(iord-1))...
                    *(kord*(REGFAC(3)^kord)*(q(3)^(kord-1)))*((REGFAC(2)^jord)*(q(2)^jord));
                end
                ddW(3,1)=ddW(1,3);
                if jord==1 && kord>1
                    ddW(2,3)=ddW(2,3)+PCOE(count)*REGFAC(2)*(kord*(REGFAC(3)^kord)...
                    *(q(3)^(kord-1)))*((REGFAC(1)^iord)*(q(1)^iord));
                elseif kord==1 && jord>1
                    ddW(2,3)=ddW(2,3)+PCOE(count)*REGFAC(3)*(jord*(REGFAC(2)^jord)...
                    *(q(2)^(jord-1)))*((REGFAC(1)^iord)*(q(1)^iord));
                elseif kord==1 && jord==1
                    ddW(2,3)=ddW(2,3)+PCOE(count)*REGFAC(3)*REGFAC(2)*((REGFAC(1)^iord)...
                    *(q(1)^iord));
                elseif kord-1<0 || jord-1<0
                    ddW(2,3)=ddW(2,3)+0;
                else
                    ddW(2,3)=ddW(2,3)+PCOE(count)*jord*(REGFAC(2)^jord)*(q(2)^(jord-1))...
                    *(kord*(REGFAC(3)^kord)*(q(3)^(kord-1)))*((REGFAC(1)^iord)*(q(1)^iord));
                end
                ddW(3,2)=ddW(2,3);
                count=count+1;
            end
            Pord(1)=Pord(1)-1;
        end
        Pord(2)=Pord(2)-1;
    end
elseif FUN==3 % reduziertes Polynom
    count=1;
    %Teilfunktion f(x,y)
    for jord=ORDER(2):-1:0
        for iord=ORDER(1):-1:0
            W=W+PCOE(count)*x(iord+1)*y(jord+1);
            dW(1)=dW(1)+PCOE(count)*dx(iord+1)*y(jord+1);
            dW(2)=dW(2)+PCOE(count)*x(iord+1)*dy(jord+1);
            ddW(1,1)=ddW(1,1)+PCOE(count)*ddx(iord+1)*y(jord+1);
            ddW(1,2)=ddW(1,2)+PCOE(count)*dx(iord+1)*dy(jord+1);
            ddW(2,1)=ddW(2,1)+PCOE(count)*dx(iord+1)*dy(jord+1);
            ddW(2,2)=ddW(2,2)+PCOE(count)*x(iord+1)*ddy(jord+1);
            count=count+1;
        end
    end
    %Teilfunktion f(x,z)
    for kord=ORDER(3):-1:1
        for iord=ORDER(1):-1:0
            W=W+PCOE(count)*x(iord+1)*z(kord+1);
            dW(1)=dW(1)+PCOE(count)*dx(iord+1)*z(kord+1);
            dW(3)=dW(3)+PCOE(count)*x(iord+1)*dz(kord+1);
            ddW(1,1)=ddW(1,1)+PCOE(count)*ddx(iord+1)*z(kord+1);
            ddW(1,3)=ddW(1,3)+PCOE(count)*dx(iord+1)*dz(kord+1);
            ddW(3,1)=ddW(3,1)+PCOE(count)*dx(iord+1)*dz(kord+1);
            ddW(3,3)=ddW(3,3)+PCOE(count)*x(iord+1)*ddz(kord+1);
            count=count+1;
        end
    end
    %Teilfunktion f(y,z)
    for kord=ORDER(3):-1:1
        for jord=ORDER(2):-1:1
            W=W+PCOE(count)*y(jord+1)*z(kord+1);
            dW(2)=dW(2)+PCOE(count)*dy(jord+1)*z(kord+1);
            dW(3)=dW(3)+PCOE(count)*y(jord+1)*dz(kord+1);
            ddW(2,2)=ddW(2,2)+PCOE(count)*ddy(jord+1)*z(kord+1);
            ddW(2,3)=ddW(2,3)+PCOE(count)*dy(jord+1)*dz(kord+1);
            ddW(3,2)=ddW(3,2)+PCOE(count)*dy(jord+1)*dz(kord+1);
            ddW(3,3)=ddW(3,3)+PCOE(count)*y(jord+1)*ddz(kord+1);
            count=count+1;
        end
    end
end
if TYPE==2                          % Berechnung der Inverse
    ddW(1,1)=((-ddW(1,1)*W^2)+(2*W*dW(1)^2))/W^4;
    ddW(1,2)=((-ddW(1,2)*W^2)+(2*W*dW(1)*dW(2)))/W^4;
    ddW(1,3)=((-ddW(1,3)*W^2)+(2*W*dW(1)*dW(3)))/W^4;
    ddW(2,2)=((-ddW(2,2)*W^2)+(2*W*dW(2)^2))/W^4;
    ddW(2,1)=((-ddW(2,1)*W^2)+(2*W*dW(2)*dW(1)))/W^4;
    ddW(2,3)=((-ddW(2,3)*W^2)+(2*W*dW(2)*dW(3)))/W^4;
    ddW(3,3)=((-ddW(3,3)*W^2)+(2*W*dW(3)^2))/W^4;
    ddW(3,1)=((-ddW(3,1)*W^2)+(2*W*dW(3)*dW(1)))/W^4;
    ddW(3,2)=((-ddW(3,2)*W^2)+(2*W*dW(3)*dW(2)))/W^4;
    dW(1)=-dW(1)/(W^2);
    dW(2)=-dW(2)/(W^2);
    dW(3)=-dW(3)/(W^2);
    W=1/W;
end
%------------------------------------ Ausgabe
F=W;
dF=dW';
ddF=ddW;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    ENDE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
