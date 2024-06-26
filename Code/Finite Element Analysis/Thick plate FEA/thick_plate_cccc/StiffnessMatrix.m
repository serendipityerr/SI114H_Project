function[K]=StiffnessMatrix(material,mesh,quadrature,ShapeOption)
%求单元刚度矩阵
K=zeros(mesh.nn*3);
%弯曲刚度矩阵
D_f=material.E/(1-material.v^2).*[1 material.v 0 
                                  material.v 1 0 
                                  0 0 (1-material.v)/2];%bending    
%剪切刚度矩阵
D_t=[material.G 0
     0 material.G];%shear
alpha=5/6;% shear stress non-uniform modification factor
beta=material.t^3/12; %integral(z^2,-t/2,t/2)
%形函数
if strcmp(ShapeOption,'Q4')
    shape_order=4;
elseif strcmp(ShapeOption,'Q8')
    shape_order=8;
elseif strcmp(ShapeOption,'Q9')
    shape_order=9;
end
for iel=1:mesh.ne
    KElem=zeros(shape_order*3);
    %高斯积分
    for ig=1:size(quadrature.points,1)
        shapeFunction=ShapeFunction(quadrature.points(ig,1),quadrature.points(ig,2),ShapeOption);%型函数在积分点处的值
        elemCoordinates=[mesh.Nid(mesh.Eid(iel,2:end),2), mesh.Nid(mesh.Eid(iel,2:end),3)];%四个节点的坐标
        jacobian=Jacobian(shapeFunction,elemCoordinates);%2*2矩阵
        Nd=jacobian.globalDerivatives;%型函数的物理坐标导数
        B_f=zeros(3,shape_order*3);%bending
        B_t=zeros(2,shape_order*3);%shear
        for ib=1:shape_order
            temp=[0 Nd(ib,1) 0         %
                  0 0 Nd(ib,2)
                  0 Nd(ib,2) Nd(ib,1)];
            B_f(:,ib*3-2:ib*3)=temp;%bending
    
            temp=[Nd(ib,1) shapeFunction.fun(ib) 0
                  Nd(ib,2) 0 shapeFunction.fun(ib)];
            B_t(:,ib*3-2:ib*3)=temp;  %shear
        end
        KElem=KElem+quadrature.weights(ig)*(beta*B_f'*D_f*B_f+alpha*material.t*B_t'*D_t*B_t)*det(jacobian.matrix);% 
    end
    ntot=mesh.Eid(iel,2:end)*3-2;%index of DOF
    %assemble to golobal stiffness matrix
    for j=1:numel(ntot)
        for k=1:numel(ntot)
            K(ntot(j):ntot(j)+2,ntot(k):ntot(k)+2) = K(ntot(j):ntot(j)+2,ntot(k):ntot(k)+2) + KElem(j*3-2:j*3,k*3-2:k*3);
        end
    end
    
end
return