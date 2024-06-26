function [d]=solveq(K,f,bc)
% a=solveq(K,f)
% [a]=solveq(K,f,bc)
%-------------------------------------------------------------
%	f : global load vector, dim(f)= nd x 1
%	bc : boundary condition matrix
%	dim(bc)= nbc x 2, nbc : number of b.c.'s
% OUTPUT:	a : solution including boundary values
%	dim(a)= nd x 1, nd : number of dof's
%-------------------------------------------------------------
[nd,nd]=size(K);
fdof=[1:nd]';
d=zeros(size(fdof));
pdof=bc(:,1);
dp=bc(:,2);
fdof(pdof)=[];
%K*u=F
s=K(fdof,fdof)\(f(fdof)-K(fdof,pdof)*dp);
d(pdof)=dp; 
d(fdof)=s;
%--------------------------end--------------------------------
