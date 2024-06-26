function [K]=assem(edof,K,Ke)
% K=assem(edof,K,Ke)
% [K]=assem(edof,K,Ke)
%-------------------------------------------------------------
% PURPOSE
%	Assemble element matrices Ke into the global
%	stiffness matrix K
%	according to the topology matrix edof.
% INPUT: edof:	dof topology matrix
%	K :	the global stiffness matrix
%	Ke:	element stiffness matrix
% OUTPUT: K :	the new global stiffness matrix
%-------------------------------------------------------------
[nie,n]=size(edof);
t=edof(:,2:n); 
for i = 1:nie
K(t(i,:),t(i,:)) = K(t(i,:),t(i,:))+Ke(:,:,i);
end
%--------------------------end--------------------------------
