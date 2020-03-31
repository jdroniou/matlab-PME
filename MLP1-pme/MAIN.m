%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass-lumped P1 code for u_t-Delta u^m=0 with Dirichlet BC
%
%    Author: Jerome Droniou
%    Date: 04/01/20
%
% This code is provided as a companion of the article
%   "The gradient discretisation method for slow and fast diffusion porous media equations", J. Droniou and K. N. Le,
%   to appear in SIAM J. Numer. Anal. https://arxiv.org/abs/1905.01785
%
%
% Usage and modification of this code is permitted, but any scientific publication resulting from
% part of this code should cite the aforementioned article.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
format long;

% Parameters for nonlinearity and barenblatt solution
global t0;
global m;
global CB;

forcediff=1;

t0=.1;
m=4;
CB=0.005;

%t0=.5;
%m=.7;
%CB=0.1;

% Final times
T=1;
dt_initial=0.1;

% nonlinear iterations
tol=1e-8;
itermax=100;
% Relaxation parameter (0 for no relaxation)
relax = 0;

%%
% Sequence of meshes over which we want to run the scheme
%%
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'mesh1_3.mat';'mesh1_4.mat'};%'mesh1_5.mat';'mesh1_6.mat'};%'mesh1_7.mat'};

nbmeshes=size(meshes,1);
Lmp1_error=zeros(nbmeshes,1);
L2beta_error=zeros(nbmeshes,1);
h=zeros(nbmeshes,1);
Ndt=zeros(nbmeshes,1);

% To see the results printed in file
fid = fopen('results.txt','w');
str = sprintf('m=%f, t0=%f\n',m,t0);
forkprint(fid,str);
%%%fclose(fid);
Ndt(1) = ceil(T/dt_initial);

for imesh=1:nbmeshes
  % Load mesh
  loadmesh=strcat('load ../../HHO-Lapl-OM/matlab_meshes/',meshes{imesh});
  str = sprintf('%s\n',loadmesh);
  forkprint(fid,str);
  eval(loadmesh);
  % Compute real centers of gravity, mesh diameter and area of dual mesh
  cg=gravity_centers(ncell,cell_v,vertex,area);
  h(imesh)=max(abs(diam));
  dualarea=compute_dualarea(area,ncell,nvert,cell_v);

  % Time steps
%  Ndt(imesh)=ceil(T/h(imesh)^2);
  if (imesh>1) 
    Ndt(imesh) = Ndt(imesh-1)*4;
  end;

  str = sprintf('mesh= %s, h= %4.2e, time step= %4.2e \n',meshes{imesh},h(imesh),T/Ndt(imesh));
  forkprint(fid,str);

  %% Initialise RHS and unknown
  %    m>1: the unknowns are the solution 'u'
  %    m<1: the unknowns are u^m
  %
  % Initial condition and exact solution
  ex_sol=exact_solution(0,vertex)';

  write_solution_vtk(ex_sol,strcat('VTKout/solution0'),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
  write_solution_vtk(ex_sol,'VTKout/ex_sol0',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

  X=transform(ex_sol,'to_unknowns');

  %% ASSEMBLE MATRIX of Laplacian
  %
  [A,b]=assemble_diffusion_system(cell_v,cell_n,ncell,nvert,vertex);

  % Time steppings
  dt=T/Ndt(imesh);
  ave_newton(imesh) = 0;
  for idt=1:Ndt(imesh);
    str = sprintf('idt=%d / %d\n',idt,Ndt(imesh));
    forkprint(fid,str);
    % Solution: non-linear iterations
    Xprev = X;

    iter = 0;
    res = 1;

    uprev = transform(Xprev,'to_solution');
    rhs = dualarea.*uprev + dt*b;

    % Dirichlet BC
    for i=1:ncell 
      I=find(cell_n{i}==0);
      if (size(I,2)>0)
        % I and I+1 are vertices on the boundary
        bdry_vert_indices = [cell_v{i}(I) cell_v{i}(I+1)];
        bdry_vert = vertex(bdry_vert_indices,:);
        rhs(bdry_vert_indices) = dualarea(bdry_vert_indices).*exact_solution(idt*dt,bdry_vert)' ...
              + forcediff*dt .*(exact_solution(idt*dt,bdry_vert).^m)';
      end
    end

    while (iter < itermax && res > tol)
      if (m>1)
        %%% Newton
        Mass = sparse(diag(dualarea)); 
        Nlin = sparse(diag(m*abs(Xprev').^(m-1)));
        Aglob = Mass +  forcediff*dt*A*Nlin;
        betau = Xprev.^m;
        nlsource = Mass * Xprev + forcediff*dt*A*betau;

        % bicgstab
        [L,U] = ilu(Aglob,struct('type','ilutp','droptol',1e-6));
        [deltaX,flag]=bicgstab(Aglob,rhs-nlsource,1e-6,20,L,U);
        if (flag ~= 0)
          flag
          error('bicgstab did not converge')
        end
%         deltaX = Aglob\(rhs-nlsource); 
        X = Xprev + (1-relax)*deltaX;
      else
        % m<1
        %%% Newton
        mu=1/m;
        Mass = sparse(diag(dualarea.*mu.*abs(Xprev).^(mu-1)));
        Aglob = Mass +  forcediff*dt*A;
        nlsource = diag(dualarea)*(abs(Xprev).^mu.*sign(Xprev)) + forcediff*dt*A*Xprev;

        % bicgstab
        [L,U] = ilu(Aglob,struct('type','ilutp','droptol',1e-6));
        [deltaX,flag]=bicgstab(Aglob,rhs-nlsource,1e-6,20,L,U);
        if (flag ~= 0)
          flag
          error('bicgstab did not converge')
        end
        X = Xprev + (1-relax)*deltaX;
          
      end

      iter = iter+1;
      % residual by increments
%      res = norm(X-Xprev,Inf) / norm(Xprev,Inf);
      % residual of non-linear system
      if (m > 1)
        betau = abs(X).^m.*sign(X);
        res = norm( (Mass*X+forcediff*dt*A*betau)-rhs , 2);
      else
        Mass = sparse(diag(dualarea));
        Xom = abs(X).^(1/m).*sign(X);
        res = norm( Mass * Xom + forcediff*dt*A*X - rhs , 2);
      end

      Xprev = X;
    end; % end nonlinear iterations
    ave_newton(imesh) = ave_newton(imesh) + iter;

   if (iter==itermax)
      res
      iter
      error('no convergence')
   end;
  usol = transform(X,'to_solution');

  % Write the solution and grid vtk files, to be plotted by "paraview"
  write_solution_vtk(usol,strcat('VTKout/solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

  % Exact solution
  ex_sol=exact_solution(dt*idt,vertex)';
  write_solution_vtk(ex_sol,strcat('VTKout/ex_sol',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
  str = sprintf('Solution computed, iter=%d, res=%4.2e, max sol=%f, max ex_sol=%f\n', iter, res, max(usol), max(ex_sol));
  forkprint(fid,str);

  end; % end time stepping

ave_newton(imesh) = ave_newton(imesh)/Ndt(imesh);

% compute error
[Lmp1_error(imesh) L2beta_error(imesh) H1beta_error(imesh)] = compute_errors(usol,ex_sol,dualarea,A);

str = sprintf('Mesh %i. Errors: L^(m+1)=%4.2e, L^2 on u^m:%4.2e, H1 on u^m:%4.2e\n',imesh,Lmp1_error(imesh),L2beta_error(imesh),H1beta_error(imesh));
forkprint(fid,str);

% fraction of negative mass
fraction_neg_mass(imesh) = abs( sum(dualarea.*min(usol,0)) / sum(dualarea.*abs(usol)) );
str = sprintf('Mesh %i. Fraction negative mass %f\n',imesh,fraction_neg_mass(imesh));
forkprint(fid,str);

% nb of interior vertices
nvert_int(imesh) = nvert;
for i=1:ncell
  for j=1:size(cell_e{i},2)
	  if (cell_n{i}(j)==0) 
      nvert_int(imesh) = nvert_int(imesh)-1;
    end
  end
end

end; % end meshes


for imesh=1:nbmeshes-1 % convergence rate
  ocLmp1(imesh)=log(Lmp1_error(imesh)/Lmp1_error(imesh+1))/log(h(imesh)/h(imesh+1));
  ocL2beta(imesh)=log(L2beta_error(imesh)/L2beta_error(imesh+1))/log(h(imesh)/h(imesh+1));
  ocH1beta(imesh)=log(H1beta_error(imesh)/H1beta_error(imesh+1))/log(h(imesh)/h(imesh+1));
end

str = sprintf('\nErrors in L^(m+1) and orders of convergence:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
  if (imesh==1)
    str = sprintf('\t%4.2e\n',Lmp1_error(imesh));
    forkprint(fid,str);
  else
    str = sprintf('\t%4.2e \t %4.2e\n',Lmp1_error(imesh),ocLmp1(imesh-1));
    forkprint(fid,str);
  end
end

str = sprintf('\nErrors in L^2 on u^m and orders of convergence:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
  if (imesh==1)
    str = sprintf('\t%4.2e\n',L2beta_error(imesh));
    forkprint(fid,str);
  else
    str = sprintf('\t%4.2e \t %4.2e\n',L2beta_error(imesh),ocL2beta(imesh-1));
    forkprint(fid,str);
  end
end

str = sprintf('\nErrors in H1 on u^m and orders of convergence:\n');
forkprint(fid,str);
for imesh=1:nbmeshes
  if (imesh==1)
    str = sprintf('\t%4.2e\n',H1beta_error(imesh));
    forkprint(fid,str);
  else
    str = sprintf('\t%4.2e \t %4.2e\n',H1beta_error(imesh),ocH1beta(imesh-1));
    forkprint(fid,str);
  end
end

fclose(fid);

% Write data file
fid = fopen('data_rates.dat','w');
fprintf(fid,'meshsize timestep L2error_zeta Lmp1error H1error NvertInt AveNewton FractionNegMass\n');
for i=1:nbmeshes
  fprintf(fid,'%f %f %f %f %f %d %f %f\n',h(i),T/Ndt(i),L2beta_error(i),Lmp1_error(i),H1beta_error(i),nvert_int(i),ave_newton(i),fraction_neg_mass(i));
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function forkprint(fid,str);

fprintf(fid,'%s',str);
fprintf(str);


end

