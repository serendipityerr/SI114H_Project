function [D]=hooke(E,v)
% D=hooke(E,v)
%-------------------------------------------------------------
% INPUT:	E : Young's modulus
%	v : Poissons const.
% OUTPUT: D : material matrix
%------------------------------------------------------------- 
D=E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];
%--------------------------end--------------------------------
