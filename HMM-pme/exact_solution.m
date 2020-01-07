% Exact solution: Barenblatt translated by t0
function u=exact_sol(t,z)

global m;
global t0;
global CB;

alpha = 2/(2*(m-1)+2);
rho = alpha/2;
gamma = alpha*(m-1)/(4*m);

sz=size(z,1);
x=z(:,1);
y=z(:,2);

normz2=(x-.5).^2+(y-.5).^2;
tps=t+t0;

for i=1:size(z,1);
  u(i)=tps.^(-alpha)*max(0, CB-gamma*normz2(i)*tps.^(-2*rho)).^(1/(m-1));
end;

