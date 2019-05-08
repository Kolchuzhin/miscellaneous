% 20.02.2015
% function test_Gt_to_Prony
%-------------------------------------------------------------------------%
t=logspace(-3,2)'; Gt=1+9*exp(-t.^0.5); 
Nd=2;
Gns=Gt_to_Prony(Gt,t,Nd);


tau=Gns(:,1); gn=Gns(:,2);
Ge=min(Gt)*0.98;
X = exp(-kron(t,1./tau'));
G_fit = Ge+X*gn;


loglog(t,G_fit,'-','LineWidth',2);
hold on;
loglog(t,Gt,'o','MarkerFaceColor','c');

legend('Prony series','data');
grid on;

xlabel('time'); 
ylabel('G(t)');
%-------------------------------------------------------------------------%
