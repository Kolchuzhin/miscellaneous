function PCOE=ROMFIT3(DATA,ORDER,FUN,TYPE)
%/////////////////////////////////////////////////////////////////////////%
% PURPOSE:
%           A.1 MATLAB-Funktion zur numerischen Ermittlung der
%           Polynomkoeffizienten
% INPUT:
%            
% OUTPUT:
%
% REFERENCE:
%           [Dr. F.Bennini, Dissertation, TU Chemnitz, 2005,
%                http://archiv.tu-chemnitz.de/pub/2005/0136]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DATA:  Matrix(:,1:4) --> DATA=[q1 q2 q3 f(q1,q2,q3)]                  %
%   ORDER: Vektor(1:3)   --> ORDER=[Nx Ny Nz]                             %
%   FUN:   ganze Zahl    --> 1 ... Vollbesetztes Polynom                  %
%                            2 ... Pascal Polynom                         %
%                            3 ... Reduziertes Polynom                    %
%   TYPE:  ganze zahl    --> 1 ... Normal                                 %
%                            2 ... Inverse                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% Tested by Kolchuzhin V.A., LMGT, TU Chemnitz, 03.09.2007 16:47
% rev. 29.01.2011 17:10 PCOE=ROMFIT3(DATA,ORDER,FUN,TYPE)
%/////////////////////////////////////////////////////////////////////////%
%=========================================================================%
if FUN==1 % vollbesetztes Polynom
    count=1;
    for zi=0:ORDER(3)
        for yi=0:ORDER(2)
            for xi=0:ORDER(1)
                PPROD(:,count)=DATA(:,1).^(ORDER(1)-xi).* ...
                DATA(:,2).^(ORDER(2)-yi).*DATA(:,3).^(ORDER(3)-zi);
                count=count+1;
            end
        end
    end
    if TYPE==1 % Normal
        PCOE=(PPROD\DATA(:,4));
    elseif TYPE==2 % Inverse
        PCOE=(PPROD\(1./DATA(:,4)));
    end
elseif FUN==2 % Pascal Polynom
%BEDINGUNG für diese schleife x>=y>=z
    count=1;
    for zi=0:ORDER(3)
        P_ORDER(1)=ORDER(1)-zi;
        P_ORDER(2)=ORDER(2);
        for yi=0:P_ORDER(2)
            for xi=0:P_ORDER(1)
                PPROD(:,count)=DATA(:,1).^(xi).*DATA(:,2).^(yi).*DATA(:,3).^(zi);
                count=count+1;
            end

    P_ORDER(1)=P_ORDER(1)-1;
        end
        P_ORDER(2)=P_ORDER(2)-1;
    end
    if TYPE==1 % Normal
        PCOE=(PPROD\DATA(:,4));
    elseif TYPE==2 % Inverse
        PCOE=(PPROD\(1./DATA(:,4)));
    end
elseif FUN==3 % reduziertes Polynom
    count=1;
    for yi=0:ORDER(2)
        for xi=0:ORDER(1)
            PRODXY(:,count)=DATA(:,1).^(ORDER(1)-xi).*DATA(:,2).^(ORDER(2)-yi);
            count=count+1;
        end
    end
    count=1;
    for zi=0:(ORDER(3)-1)
        for xi=0:ORDER(1)
            PRODXZ(:,count)=DATA(:,1).^(ORDER(1)-xi).*DATA(:,3).^(ORDER(3)-zi);
            count=count+1;
        end
    end
    count=1;
    for zi=0:(ORDER(3)-1)
        for yi=0:(ORDER(2)-1)
            PRODYZ(:,count)=DATA(:,2).^(ORDER(2)-yi).*DATA(:,3).^(ORDER(3)-zi);
            count=count+1;
        end
    end
    PPROD=[PRODXY PRODXZ PRODYZ];
    if TYPE==1 % Normal
        PCOE=(PPROD\DATA(:,4));
    elseif TYPE==2 % Inverse
        PCOE=(PPROD\(1./DATA(:,4)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     ENDE     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
