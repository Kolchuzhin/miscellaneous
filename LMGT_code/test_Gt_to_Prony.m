% 20.02.2015
%-------------------------------------------------------------------------%
function test_Gt_to_Prony

t=logspace(-3,2)'; Gt=1+9*exp(-t.^0.5); 
Nd=2;
Gns=Gt_to_Prony(Gt,t,Nd);

loglog(t,Gt,'-','LineWidth',2);


tau=Gns(:,1); gn=Gns(:,2);
Ge=min(Gt)*0.98;
X = exp(-kron(t,1./tau'));
G_fit = Ge+X*gn;

hold on;
loglog(t,G_fit,'o','MarkerFaceColor','c');
legend('data','Prony series');
grid on;

xlabel('time'); 
ylabel('G(t)');
%-------------------------------------------------------------------------%
