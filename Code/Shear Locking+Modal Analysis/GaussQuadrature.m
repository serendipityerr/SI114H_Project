
function[quadrature]=GaussQuadrature(GaussOption)
% Compute the Gauss Point with their weights
% Case 1: Quad 1x1
% Case 2: Quad 2x2
% Case 3: Quad 3x3
% Gauss points are automatically chosen relative to shape function

quadrature=struct();

switch GaussOption
    case 'gauss1'
        quadrature.points=[0 0];
        
        quadrature.weights=4;
        
    case 'gauss2'
        quadrature.points=[...
            -1/sqrt(3) -1/sqrt(3)
             1/sqrt(3) -1/sqrt(3)
             1/sqrt(3)  1/sqrt(3)
            -1/sqrt(3)  1/sqrt(3)];
        
        quadrature.weights=ones(size(quadrature.points,1),1);
        
    case 'gauss3'
        quadrature.points=[...
                    0           0
             sqrt(3/5)  sqrt(3/5)
            -sqrt(3/5)  sqrt(3/5)
            -sqrt(3/5) -sqrt(3/5)
             sqrt(3/5) -sqrt(3/5)
                     0  sqrt(3/5)
            -sqrt(3/5)         0
                     0 -sqrt(3/5)
             sqrt(3/5)         0];
         
         quadrature.weights=[...
             64/81
             25/81
             25/81
             25/81
             25/81
             40/81
             40/81
             40/81
             40/81];
end

return