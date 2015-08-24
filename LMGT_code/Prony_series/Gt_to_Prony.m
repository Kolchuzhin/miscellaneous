% Simple MatLab program for Prony series approximation of G(t)
% KMB Jansen version 3-11-2004

function Gns=Gt_to_Prony(Gtt,tt,Nd)
% For Prony parameters determination of a data file
% Nd = desired logarithmic spacing of tau (points/decade) [1..5], default=2
% example: t=logspace(-3,2)'; G=1+9*exp(-t.^0.5); Gt_to_Prony(G,t,2);

% PRONY FITTING
Ge=min(Gtt)*0.98;                           % rubbery plateau
logt0=log10(min(tt)); logt1=log10(max(tt));
tau = 10.^(logt0:1/Nd:logt1)';              % chosen Prony times
% xx = kron(ww,tau'); X = xx.^2./(1+xx.^2); %   dynamic data (Gw,w)
xx = kron(tt,1./tau'); X = exp(-xx);        % transient data (Gt,t)
%gn = (X\(Gtt-Ge));         % solution by matrix inversion (same as lsqlin)
gn = lsqnonneg(X,Gtt-Ge);   % solution by non neg least squares (in MatLab 5: use nnls)
G_fit = Ge+X*gn;            % for plotting fitted values: loglog(tt,G_fit);

% OUTPUT TO FUNCTION
Gns=[tau gn];
