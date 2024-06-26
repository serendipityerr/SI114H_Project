clear;clc

%---------------------GEOMETRY-----------------------
a=1;  %length in x
b=1;  %length in y
%---------------------MESH---------------------------
ex=32; %unitx 
ey=32; %unity
mesh=MeshGenerator(a,b,ex,ey);
%---------------------ELEMENT------------------------
ShapeOption='Q4';%Q8 Q9
quadrature=GaussQuadrature('gauss2');% Q4 element points and weights
%---------------------MATERIAL------------------------
%SI mm t MPa N
material=struct();
material.t=1;%thickness
material.E=2.1*10^11;%2.1e5MPa
material.v=0.3;%Poison ratio
material.rho=2.7e-9;%密度
material.G=material.E/(2+2*material.v);%剪切模量

K=StiffnessMatrix(material,mesh,quadrature,ShapeOption);%K
M=MassMatrix(material,mesh,quadrature,ShapeOption);%质量矩阵
num_points = (ex + 1) * (ey + 1);
label = zeros(num_points, 4); % Initialize label matrix

for i = 1:num_points
    label(i, 1) = i; % Point number
    label(i, 2) = 3 * (i - 1) + 1; % Degree of freedom 1
    label(i, 3) = 3 * (i - 1) + 2; % Degree of freedom 2
    label(i, 4) = 3 * (i - 1) + 3; % Degree of freedom 3
end
%---------------------LOADS----------------------------
p0=1e6;
F=ForceVector(mesh,quadrature,ShapeOption,p0);%F
% F=ForceVector(cl,c0,mesh)
%---------------------1e-3-----------------------
%nc is nodes of constrained elemment
nc=[mesh.lato1 mesh.lato2(2:end) mesh.lato3(2:end) mesh.lato4(2:end-1)];
[K_c,M_c,F_c,nctot]=Constraints(nc,K,M,F);
%---------------------SOLUTION---------------------------
[w,thetax,thethy]=StaticSolver(K_c,F_c,mesh,nctot);
%---------------------PLOT-------------------------------
% figure
% surf(mesh.x(1,:),mesh.y(:,1),reshape(w,size(mesh.x,2),size(mesh.x,1))');
% % axis equal
% hold on
% title('Deformed geometry')
% %3.scatter

figure;
scatter3(mesh.x,mesh.y,reshape(w,size(mesh.x,2),size(mesh.x,1)), '.');zlim([0,0.0000025])
title('Final Displacement of Each Point of Mindlin');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Displacement (mm)');
% colorbar;
% 
% % 3.scatter
% figure;
% scatter3(mesh.x,mesh.y,reshape(w,size(mesh.x,2),size(mesh.x,1))', '*');
% title('Initial Position of Each Point');
% xlabel('X Position (mm)');
% ylabel('Y Position (mm)');
% zlabel('Initial Position');
% colorbar;
% 
% 
% figure;
% imagesc(mesh.x(1,:),mesh.y(:,1),reshape(w,size(mesh.x,2),size(mesh.x,1))');
% title('Deflection Heatmap');
% xlabel('X Position (mm)');
% ylabel('Y Position (mm)');
% colorbar;
% axis equal;
% 
% % Compute internal forces
% num_elements = size(mesh.Eid, 1);
% Shear = zeros(mesh.nn, 1);
% Mx = zeros(mesh.nn, 1);
% My = zeros(mesh.nn, 1);
% 
% for i = 1:num_elements
%     element_nodes = mesh.Eid(i, 2:end);
%     element_displacement = [w(element_nodes); thetax(element_nodes); thethy(element_nodes)];
%     [shear, mx, my] = ElementInternalForces(material, mesh, element_nodes, element_displacement);
% 
%     Shear(element_nodes) = Shear(element_nodes) + shear;
%     Mx(element_nodes) = Mx(element_nodes) + mx;
%     My(element_nodes) = My(element_nodes) + my;
% end
% 
% % Plot internal forces with patch
% % 1. Shear plot
% figure;
% title('Shear');
% patch('Faces', mesh.Eid(:,2:5), 'Vertices', [mesh.cordinates(:,1), mesh.cordinates(:,2)], 'FaceVertexCData', Shear, 'EdgeColor', 'none', 'FaceColor', 'interp');
% colorbar;
% 
% % 2. Mx plot
% figure;
% title('Mx');
% patch('Faces', mesh.Eid(:,2:5), 'Vertices', [mesh.cordinates(:,1), mesh.cordinates(:,2)], 'FaceVertexCData', Mx, 'EdgeColor', 'none', 'FaceColor', 'interp');
% colorbar;
% 
% % 3. My plot
% figure;
% title('My');
% patch('Faces', mesh.Eid(:,2:5), 'Vertices', [mesh.cordinates(:,1), mesh.cordinates(:,2)], 'FaceVertexCData', My, 'EdgeColor', 'none', 'FaceColor', 'interp');
% colorbar;
% 
% 
% %visualize
% x=mesh.x(1,:);
% y=mesh.y(:,1);
% z=reshape(w,size(mesh.x,2),size(mesh.x,1));
% figure;
% h = surf(x, y, zeros(size(z)));
% title('Dynamic Displacement Process');
% xlabel('X Position (mm)');
% ylabel('Y Position (mm)');
% zlabel('Displacement (mm)');
% colorbar;
% zlim([0 max(z(:))]);
% % update
% num_steps = 100;
% for step = 1:num_steps
%     current_z = (step / num_steps) * z;
%     set(h, 'ZData', current_z);
%     drawnow;
%     pause(0.1); % gap
% end
% NodeTable = table((1:mesh.nn)', w, thetax, thethy, ...
%     'VariableNames', {'Node', 'Deflection', 'ThetaX', 'ThetaY'});
% % Display the table in a figure using uitable
% figure;
% uitable('Data', NodeTable{:,:}, 'ColumnName', NodeTable.Properties.VariableNames, ...
%     'Units', 'normalized', 'Position', [0 0 1 1]);
% density = 20; % 调整密度参数以减少箭头数量
% [xq, yq] = meshgrid(linspace(min(mesh.x(:)), max(mesh.x(:)), density), ...
%                     linspace(min(mesh.y(:)), max(mesh.y(:)), density));
% % 插值位移和旋转以匹配箭头位置
% uq = interp2(mesh.x, mesh.y, reshape(thetax, size(mesh.x,2), size(mesh.x,1))', xq, yq);
% vq = interp2(mesh.x, mesh.y, reshape(thethy, size(mesh.x,2), size(mesh.x,1))', xq, yq);
% % 绘制背景位移颜色渐变
% figure;
% contourf(mesh.x, mesh.y, reshape(w, size(mesh.x,2), size(mesh.x,1))', 20, 'LineColor', 'none');
% colormap jet; 
% colorbar;
% hold on;
% % 绘制矢量场
% quiver(xq, yq, uq, vq, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5); % 调整箭头颜色和大小
% title('Displacement Vector Field');
% xlabel('X Position (mm)');
% ylabel('Y Position (mm)');
% axis equal tight;
% grid on;
% scale_factor = 1; % 调整比例因子以控制箭头的长度
% quiver(xq, yq, uq*scale_factor, vq*scale_factor, 0, 'Color', 'k');
% legend({'Background Deflection', 'Displacement Vectors'}, 'Location', 'best');
% hold off;


