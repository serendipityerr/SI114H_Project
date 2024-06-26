function[shapeFunction]=ShapeFunction(xi,eta,ShapeOption)
% Database of shape functions and their natural derivatives (xi,eta)
% 'quad1': 4 nodes
% 'quad2': 8 nodes
% 'quad3': 9 nodes

shapeFunction=struct();

switch ShapeOption
    
    case 'Q4'
        shapeFunction.fun=1/4*[(1-xi)*(1-eta)
                               (1+xi)*(1-eta)
                               (1+xi)*(1+eta)
                               (1-xi)*(1+eta)];
                           
        shapeFunction.dfun=1/4*[-(1-eta),-(1-xi)
                                1-eta,-(1+xi)
                                1+eta,1+xi
                                -(1+eta),1-xi];
                            
    case 'Q8'
        shapeFunction.fun=[1/4*xi*(1-xi)*eta*(1-eta)
                           -1/2*xi*(1-xi)*(1+eta)*(1-eta)
                           -1/4*xi*(1-xi)*eta*(1+eta)
                           1/2*(1+xi)*(1-xi)*(1+eta)*eta
                           1/4*xi*(1+xi)*eta*(1+eta)
                           1/2*xi*(1+xi)*(1+eta)*(1-eta)
                           -1/4*xi*(1+xi)*eta*(1-eta)
                           -1/2*(1+xi)*(1-xi)*(1-eta)*eta];
                       
       shapeFunction.dfun=[...
           1/4*eta*(-1+eta)*(-1+2*xi)      1/4*xi*(-1+xi)*(-1+2*eta)
      -1/2*(1+eta)*(-1+eta)*(-1+2*xi)                -xi*(-1+xi)*eta
            1/4*eta*(1+eta)*(-1+2*xi)       1/4*xi*(-1+xi)*(1+2*eta)
                      -xi*eta*(1+eta)  -1/2*(1+xi)*(-1+xi)*(1+2*eta)
             1/4*eta*(1+eta)*(1+2*xi)        1/4*xi*(1+xi)*(1+2*eta)
       -1/2*(1+eta)*(-1+eta)*(1+2*xi)                 -xi*(1+xi)*eta
            1/4*eta*(-1+eta)*(1+2*xi)       1/4*xi*(1+xi)*(-1+2*eta)
                     -xi*eta*(-1+eta) -1/2*(1+xi)*(-1+xi)*(-1+2*eta)];
                       
end

return 