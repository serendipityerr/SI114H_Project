function[mesh]=MeshGenerator(a,b,ex,ey)
% Define grid
% ...

mesh=struct();

nx=ex+1; ny=ey+1;
mesh.ne=ex*ey; %number of elements
mesh.nn=nx*ny; %number of nodes

dx=a/ex; dy=b/ey;

[mesh.x,mesh.y]=meshgrid(0:dx:a,0:dy:b);
mesh.cordinates(:,1)=reshape(mesh.x',numel(mesh.x),1);%x coordinates
mesh.cordinates(:,2)=reshape(mesh.y',numel(mesh.y),1);%y coordinates

% side 1
mesh.lato1=1:nx;
% side 2
mesh.lato2=nx:nx:mesh.nn;
% side 3
mesh.lato3=mesh.nn:-1:mesh.nn-ex;
% side 4
mesh.lato4=mesh.nn-ex:-nx:1;

mesh.Nid=zeros(mesh.nn,4);% [node_id x y z]
mesh.Eid=zeros(mesh.ne,5);% [ele_id node1 node2 node3 node4]

for j=1:ny
    for i=1:nx
        mesh.Nid(i+nx*(j-1),:)=[i+nx*(j-1), (i-1)*dx, (j-1)*dy, 0];
    end
end

for j=1:ey
    for i=1:ex
        mesh.Eid(i+ex*(j-1),:)=[i+ex*(j-1), i+nx*(j-1)...
                           i+nx*(j-1)+1 i+nx*(j-1)+nx+1 i+nx*(j-1)+nx];
    end
end
return