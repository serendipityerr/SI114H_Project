function[F]=ForceVector(mesh,quadrature,ShapeOption,p0)
%....


F=zeros(mesh.nn*3,1);

for iel=1:mesh.ne
    FElem=zeros(12,1);
    for ig=1:size(quadrature.points,1)
        shapeFunction=ShapeFunction(quadrature.points(ig,1),quadrature.points(ig,2),ShapeOption);
        elemCoordinates=[mesh.Nid(mesh.Eid(iel,2:end),2), mesh.Nid(mesh.Eid(iel,2:end),3)];
        jacobian=Jacobian(shapeFunction,elemCoordinates);

        for ib=1:4
            temp=[shapeFunction.fun(ib) 0 0];
            N(:,ib*3-2:ib*3)=temp;
        end

        FElem=FElem+quadrature.weights(ig)*(N'*p0)*det(jacobian.matrix);
    end

    ntot=mesh.Eid(iel,2:end)*3-2;%DOF index of whole model

    for j=1:numel(ntot)
        F(ntot(j):ntot(j)+2)=F(ntot(j):ntot(j)+2)+FElem(j*3-2:j*3);
    end 
end


return