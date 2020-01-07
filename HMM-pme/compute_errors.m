% Compute the L2 error on the function and the gradient
%
function [Lmp1,L2beta]=compute_errors(Y,Yexact,weights)

global m;

norm_mp1_exact = sum(weights.*abs(Yexact).^(m+1))^(1/(m+1));
error_mp1 = sum(weights.*abs(Yexact-Y).^(m+1))^(1/(m+1));
Lmp1 = error_mp1 / norm_mp1_exact;

norm_L2beta_exact = sum(weights.*abs(Yexact).^(2*m))^(1/2);
sYexact=sign(Yexact);
sY=sign(Y);
error_L2beta = sum(weights.*(abs(Yexact).^m.*sYexact-abs(Y).^m.*sY).^2 )^(1/2);
L2beta = error_L2beta / norm_L2beta_exact;



