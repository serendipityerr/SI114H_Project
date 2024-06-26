function [ed]=extract(edof,a)
% ed=extract(edof,a)
%-------------------------------------------------------------
% INPUT:	a:	the global displacement vector
%	edof:	topology matrix
% OUTPUT: ed:	element displacement matrix
%-------------------------------------------------------------
[nie,n]=size(edof);
%
t=edof(:,2:n);
%
for i = 1:nie
ed(i,1:(n-1))=a(t(i,:))';
end


%--------------------------end--------------------------------
