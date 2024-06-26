function [Edof]=caldof(Enode,ndof)
% Edof=caldof(Enode,ndof)
%-------------------------------------------------------------
% INPUT:	Enode:	element node
%	ndof :	degrees of freedom in each node
% OUTPUT: Edof : topology of the structure
%-------------------------------------------------------------
Nele=size(Enode,1); 
Esize=size(Enode,2)-1; 
for i=1:Nele
  Edof(i,1)=i; 
    for j=1:Esize
    Edof(i,3*j-1)=3*(Enode(i,j+1)-1)+1;
    Edof(i,3*j)=3*(Enode(i,j+1)-1)+2; Edof(i,3*j+1)=3*(Enode(i,j+1)-1)+3;
    end
end
%--------------------------end--------------------------------
