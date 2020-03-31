% Compute the L2 error on the function u and the gradient of u^m
%
function [Lmp1,L2beta,H1beta]=compute_errors(ncell,Y,Yexact,weights,DiffMat)

global m;

Y_cell = Y(1:ncell);
Yexact_cell = Yexact(1:ncell);

norm_mp1_exact = sum(weights.*abs(Yexact_cell).^(m+1))^(1/(m+1));
error_mp1 = sum(weights.*abs(Yexact_cell-Y_cell).^(m+1))^(1/(m+1));
Lmp1 = error_mp1 / norm_mp1_exact;

norm_L2beta_exact = sum(weights.*abs(Yexact_cell).^(2*m))^(1/2);
sYexact=sign(Yexact_cell);
sY=sign(Y_cell);
error_L2beta = sum(weights.*(abs(Yexact_cell).^m.*sYexact-abs(Y_cell).^m.*sY).^2 )^(1/2);
L2beta = error_L2beta / norm_L2beta_exact;

% Energy
betaY = sign(Y) .* abs(Y).^m; 
betaYexact = sign(Yexact) .* abs(Yexact).^m;
error_H1beta = ((betaY-betaYexact)' * DiffMat * (betaY-betaYexact) )^(1/2);
norm_H1beta_exact = ((betaYexact)' * DiffMat * (betaYexact) )^(1/2);
H1beta = error_H1beta / norm_H1beta_exact;

