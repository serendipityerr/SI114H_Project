
function[M]=MassMatrix(material,mesh,quadrature,ShapeOption)

M=zeros(mesh.nn*3);

I=[material.t 0 0 %w
   0 material.t^3/12 0 %dw/dx 转动惯量*角加速度
   0 0 material.t^3/12];%dw/dy
%形函数类型定义
if strcmp(ShapeOption,'Q4')
    shape_order=4;
elseif strcmp(ShapeOption,'Q8')
    shape_order=8;
elseif strcmp(ShapeOption,'Q9')
    shape_order=9;
end

for iel=1:mesh.ne
    MElem=zeros(shape_order*3);
    %高斯积分
    for ig=1:size(quadrature.points,1)
        shapeFunction=ShapeFunction(quadrature.points(ig,1),quadrature.points(ig,2),ShapeOption);
        elemCoordinates=[mesh.Nid(mesh.Eid(iel,2:end),2), mesh.Nid(mesh.Eid(iel,2:end),3)];
        jacobian=Jacobian(shapeFunction,elemCoordinates);
        
        N=zeros(3,shape_order*3);

        for ib=1:shape_order
            temp=[shapeFunction.fun(ib) 0 0
                  0 shapeFunction.fun(ib) 0
                  0 0 shapeFunction.fun(ib)];
            N(:,ib*3-2:ib*3)=temp;
        end
        MElem=MElem+quadrature.weights(ig)*(material.rho*N'*I*N)*det(jacobian.matrix);
    end
    %将单元质量矩阵添加到全质量矩阵
    ntot=mesh.Eid(iel,2:end)*3-2;
    for j=1:numel(ntot)
        for k=1:numel(ntot)
            M(ntot(j):ntot(j)+2,ntot(k):ntot(k)+2) = M(ntot(j):ntot(j)+2,ntot(k):ntot(k)+2) + MElem(j*3-2:j*3,k*3-2:k*3);
        end
    end
    
end
return
