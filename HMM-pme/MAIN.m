%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HMM code for u_t-Delta u^m=0 with Dirichlet BC
%
%		Author: Jerome Droniou
%		Date: 04/01/20
%
% This code is provided as a companion of the article
%   "The gradient discretisation method for slow and fast diffusion porous media equations", J. Droniou and K. N. Le,
%   https://arxiv.org/abs/1905.01785
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
m=3;
CB=0.005;

%t0=.5;
%m=.3;
%CB=0.1;

% Final times
T=1;

% nonlinear iterations
tol=1e-8;
itermax=100;
% Relaxation parameter (0 for no relaxation)
relax = 0;

%%
% Sequence of meshes over which we want to run the scheme
%%
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'hexa1_1.mat';'hexa1_2.mat'};%'hexa1_3.mat';'hexa1_4.mat'};

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

for imesh=1:nbmeshes
	% Load mesh
	loadmesh=strcat('load ../../HHO-Lapl-OM/matlab_meshes/',meshes{imesh});
	str = sprintf('%s\n',loadmesh);
  forkprint(fid,str);
	eval(loadmesh);
	% Compute real centers of gravity and mesh diameter
	cg=gravity_centers(ncell,cell_v,vertex,area);
  h(imesh)=max(abs(diam));

  % Time steps
  Ndt(imesh)=ceil(T/h(imesh)^2);

  str = sprintf('mesh= %s, h= %4.2e, time step= %4.2e \n',meshes{imesh},h(imesh),T/Ndt(imesh));
  forkprint(fid,str);

	%% Initialise RHS and unknown
	%		m>1: the unknowns are the solution 'u' in the cells (1:ncell) and 'u^m' on the edges (ncell+1:ncell+nedge)
  %   m<1: the unknowns are u^m everywhere
  ex_sol=zeros(ncell+nedge,1);
	ex_sol(1:ncell)=exact_solution(0,cg);
  for i=1:ncell
		nbe=size(cell_e{i},2);
		% midpoints of edges
		vertex_loc=vertex(cell_v{i},:);
		xs=zeros(nbe,2);
		xs([1:nbe],:)=(vertex_loc([1:nbe],:)+vertex_loc([2:nbe+1],:))/2;
		ex_sol(ncell+cell_e{i})=exact_solution(0,xs);
	end;

  write_solution_vtk(ex_sol,strcat('VTKout/solution0'),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
  write_solution_vtk(ex_sol,'VTKout/ex_sol0',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

  X=transform(ex_sol,ncell,'to_unknowns');

	%% ASSEMBLE MATRIX of Laplacian
	%
	[A,b]=assemble_diffusion_system(cell_v,cell_n,cell_e,ncell,nedge,vertex,area,center,cg);

  % Time steppings
  dt=T/Ndt(imesh);
  for idt=1:Ndt(imesh);
	  str = sprintf('idt=%d / %d\n',idt,Ndt(imesh));
    forkprint(fid,str);
	  % Solution: non-linear iterations
    Xprev = X;
    iter = 0;
    res = 1;

    uprev = transform(Xprev,ncell,'to_solution');
 	  rhs = [area.*uprev(1:ncell);zeros(nedge,1)] + dt*b;

    % Dirichlet BC
    %  These BC are set up based on the fact that, on rows corresponding to boundary edges,
    %   A has 1 on the diagonal and 0 elsewhere. Since it's multiplied by forcediff * dt in the
    %   global system, so is the exact solution to define the BC
  	for i=1:ncell 
		  I=find(cell_n{i}==0);
		  if (size(I,2)>0)
			  % midpoints of the boundary edges
			  xs=zeros(size(I,2),2);
			  xs=(vertex(cell_v{i}(I),:)+vertex(cell_v{i}(I+1),:))/2;
  		  rhs(ncell+cell_e{i}(I)) = forcediff*dt*(exact_solution(idt*dt,xs).^m)';
		  end;
    end

    while (iter < itermax && res > tol)
      if (m>1)
        %%% Newton
        Mass = sparse(diag([area' zeros(1,nedge)])); 
        Nlin = sparse(diag([m*abs(Xprev(1:ncell)').^(m-1) ones(1,nedge)]));
        Aglob = Mass +  forcediff*dt*A*Nlin;
        betau = Xprev;
        betau(1:ncell) = X(1:ncell).^m;
        nlsource = Mass * [Xprev(1:ncell); zeros(nedge,1)] + forcediff*dt*A*betau;

        % bicgstab
        [L,U] = ilu(Aglob,struct('type','ilutp','droptol',1e-6));
        [deltaX,flag]=bicgstab(Aglob,rhs-nlsource,1e-6,20,L,U);
        if (flag ~= 0)
          flag
          error('bicgstab did not converge')
        end
        X = Xprev + (1-relax)*deltaX;

      else
        % m<1
        %%% Newton
        mu=1/m;
        Mass = sparse(diag([area'*mu.*abs(Xprev(1:ncell)').^(mu-1) zeros(1,nedge)]));
        Aglob = Mass +  forcediff*dt*A;
        nlsource = [ diag(area)*(abs(Xprev(1:ncell)).^mu.*sign(Xprev(1:ncell))); zeros(nedge,1)] ...
                  + forcediff*dt*A*Xprev;

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
        betau = X;
        betau(1:ncell) = X(1:ncell).^m;
        res = norm( (Mass*X+forcediff*dt*A*betau)-rhs , Inf);
      else
        Mass = sparse(diag([area' zeros(1,nedge)]));
        Xom = abs(X).^(1/m).*sign(X);
        res = norm( Mass * Xom + forcediff*dt*A*X - rhs , Inf);
      end

      Xprev = X;
    end; % end nonlinear iterations

   if (iter==itermax)
      res
      iter
      error('no convergence')
   end;
  usol = transform(X,ncell,'to_solution');

  % Write the solution and grid vtk files, to be plotted by "paraview"
  write_solution_vtk(usol,strcat('VTKout/solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

  % Exact solution
  ex_sol=zeros(ncell+nedge,1);
  ex_sol(1:ncell)=exact_solution(dt*idt,cg);
	  for i=1:ncell
		  nbe=size(cell_e{i},2);
		  % midpoints of edges
		  vertex_loc=vertex(cell_v{i},:);
		  xs=zeros(nbe,2);
		  xs([1:nbe],:)=(vertex_loc([1:nbe],:)+vertex_loc([2:nbe+1],:))/2;
		  ex_sol(ncell+cell_e{i})=exact_solution(dt*idt,xs);
	  end;
  write_solution_vtk(ex_sol,strcat('VTKout/ex_sol',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
  str = sprintf('Solution computed, iter=%d, res=%4.2e, max sol=%f, max ex_sol=%f\n', iter, res, max(usol), max(ex_sol));
  forkprint(fid,str);

  end; % end time stepping

% compute error
[Lmp1_error(imesh) L2beta_error(imesh)] = compute_errors(usol(1:ncell),ex_sol(1:ncell),area);

str = sprintf('Mesh %i. Errors: L^(m+1)=%4.2e, L^2 on u^m:%4.2e\n',imesh,Lmp1_error(imesh),L2beta_error(imesh));
forkprint(fid,str);

end; % end meshes


for imesh=1:nbmeshes-1 % convergence rate
  ocLmp1(imesh)=log(Lmp1_error(imesh)/Lmp1_error(imesh+1))/log(h(imesh)/h(imesh+1));
  ocL2beta(imesh)=log(L2beta_error(imesh)/L2beta_error(imesh+1))/log(h(imesh)/h(imesh+1));
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

fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function forkprint(fid,str);

fprintf(fid,'%s',str);
fprintf(str);

end

