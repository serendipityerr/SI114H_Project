function[jacobian]=Jacobian(shapeFunction,elemCoordinates)
% Compute Jacobian and calculate global cordinates
jacobian=struct();
jacobian.matrix=elemCoordinates'*shapeFunction.dfun;
jacobian.invmatrix=inv(jacobian.matrix);
jacobian.globalDerivatives=(jacobian.invmatrix*shapeFunction.dfun')';%对物理坐标求偏导

return