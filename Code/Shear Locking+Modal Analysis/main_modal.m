clear;close all;clc
%---------------------GEOMETRY-----------------------
a=1;
b=1;
%---------------------MESH---------------------------
ex=64;
ey=64;
mesh=MeshGenerator(a,b,ex,ey);
%---------------------ELEMENT------------------------
ShapeOption='Q4';%Q8 Q9
quadrature=GaussQuadrature('gauss1');
% gauss1
%---------------------MATERIAL------------------------
%SI mm t MPa N
material=struct();
material.t=0.01;%thickness
material.E=2.1*10^11;%2.1e5MPa
material.v=0.3;%Poison ratio
material.rho=2.7e-9;
material.G=material.E/(2+2*material.v);

K=StiffnessMatrix(material,mesh,quadrature,ShapeOption);
M=MassMatrix(material,mesh,quadrature,ShapeOption);
%---------------------LOADS----------------------------
p0=1e6;
F=ForceVector(mesh,quadrature,ShapeOption,p0);
% c0=input('\n集中荷载的大小: ');
% cl=input('\nLato da caricare: ');
% F=ForceVector(cl,c0,mesh);Lato da caricare
%---------------------CONSTRAINTS-----------------------
%nc is nodes of constrained elemment
nc=[mesh.lato1 mesh.lato2(2:end) mesh.lato3(2:end) mesh.lato4(2:end-1)];%%四边固定
[K_c,M_c,F_c,nctot]=Constraints(nc,K,M,F);
%---------------------SOLUTION---------------------------
[w,thetax,thethy]=StaticSolver(K_c,F_c,mesh,nctot);
sprintf('最大位移为：%f',max(abs(w)))
numberMode=6;
[V,D]=eigs(K_c\M_c,numberMode,'lm');
scal=0.5;
f= power(diag(sqrt(D)),-1)/(2*pi); 
Vf=sortrows([V',f],size(V,1)+1);
V=Vf(:,1:size(V,1))';
index=1:mesh.nn*3;
Vo=zeros(mesh.nn*3,numberMode);
index(nctot)=[];
for i=1:numel(index)
    Vo(index(i),:)=V(i,:);
end
f=sort(f);%排序后的频率
sprintf('频率为%f',f)
figure;
for i=1:numberMode
    subplot(3,2,i);
    surf(mesh.x(1,:),mesh.y(:,1),reshape(Vo(1:3:end,i),size(mesh.x,2),size(mesh.x,1))');
    title(sprintf('modal %d',i))
    hold on
end
for i = 1:numberMode
    figure; 
    subplot(2, 1, 1);
    plot(mesh.x(1,:), reshape(Vo(1:3:end, i), size(mesh.x, 2), size(mesh.x, 1))');
    title(sprintf('Mode %d: X-Displacement', i));
    xlabel('X'); ylabel('Displacement');
    subplot(2, 1, 2);
    plot(mesh.y(:,1), reshape(Vo(2:3:end, i), size(mesh.y, 1), size(mesh.y, 2))');
    title(sprintf('Mode %d: Y-Displacement', i));
    xlabel('Y'); ylabel('Displacement');
end
figure;
for i = 1:numberMode
    subplot(3,2,i);
    contourf(mesh.x(1,:), mesh.y(:,1), reshape(Vo(1:3:end, i), size(mesh.x, 2), size(mesh.x, 1))', 'LineStyle', 'none');
    title(sprintf('Mode %d: Contour Plot of X-Displacement', i));
    colorbar;
    colormap(jet); 
    shading interp; 
end
figure;
for i = 1:numberMode
    subplot(3,2,i);
    quiver(mesh.x(1,:), mesh.y(:,1), reshape(Vo(1:3:end, i), size(mesh.x, 2), size(mesh.x, 1))', reshape(Vo(2:3:end, i), size(mesh.y, 1), size(mesh.y, 2))');
    title(sprintf('Mode %d: Displacement Vector Field', i));
end
