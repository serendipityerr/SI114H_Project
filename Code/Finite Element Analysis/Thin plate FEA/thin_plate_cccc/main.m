clear all;clc;close all;
%----- Topology -------------------------------------------------
%square plate with n element each side
n_length=64;%number of element in long side
m_width=64;%number of element in wide side
uni=1/64;%element size
q=1e+6;%1N/mm2
%q = (pi^6 * (uni)^4 / 4);
num = 0;
for j = 1:m_width
    for i = 1:n_length
        num = num + 1;
        top_left_node = (j-1)*(n_length+1) + i;
        top_right_node = top_left_node + 1;
        bottom_right_node = j*(n_length+1) + i + 1;
        bottom_left_node = bottom_right_node - 1;
        Enode(num,:) = [num, top_left_node, top_right_node, bottom_right_node, bottom_left_node];
        ex(num,:) = [(i-1)*uni, i*uni, i*uni, (i-1)*uni];
        ey(num,:) = [(j-1)*uni, (j-1)*uni, j*uni, j*uni];
    end
end
ndof=3;
Edof=caldof(Enode,ndof);%1 12
%----- Generate Label -------------------------------------------
num_points = (n_length + 1) * (m_width + 1);
label = zeros(num_points, 4); % Initialize label matrix

for i = 1:num_points
    label(i, 1) = i; % Point number
    label(i, 2) = 3 * (i - 1) + 1; % Degree of freedom 1
    label(i, 3) = 3 * (i - 1) + 2; % Degree of freedom 2
    label(i, 4) = 3 * (i - 1) + 3; % Degree of freedom 3
end
%----- Element stiffness ----------------------------------------
E=2.1*10^11; %2.1e5MPa
v=0.3;%Poison ratio
D=hooke(E,v);%stress=D*strain %% Mechanics of Plate and Shell 
t=0.01; %thickness of the plate
ep=[t];% thickness of element
nie=size(Enode,1); %number of element
for i=1:nie
Ke(:,:,i)=platre(ex(i,:),ey(i,:),ep,D); %%% Mechanics of Plate and Shell 
end

%----- Assemble Ke into K ---------------------------------------
K=zeros(max(max(Edof))); K=assem(Edof,K,Ke);

%----- load vector f and boundary conditions bc -----------------
f=zeros(max(max(Edof)),1);	
f(1:3:end)=-q*uni*uni;%
num=0;
%find boundary node
for i=0:m_width 
    num=num+1;
    boundary_node(num)=(i)*(n_length+1)+1;
    num=num+1;
    boundary_node(num)=(i+1)*(n_length+1);
end
boundary_node(num+1:num+1+n_length)=1:n_length+1;
boundary_node(num+2+n_length:num+2+2*n_length)=(m_width)*(n_length+1)+1:(m_width+1)*(n_length+1);
%map boundary constraint to edof(element degree of freedom)
for i=1:length(boundary_node)
    bc((i-1)*3+1,:)=[3*(boundary_node(i)-1)+1 0]; %0 is displacement
    bc((i-1)*3+2,:)=[3*(boundary_node(i)-1)+2 0];
    bc((i-1)*3+3,:)=[3*(boundary_node(i)-1)+3 0];
%     bc(i,:)=[3*(boundary_node(i)-1)+1 0]; %simple supported
end
%----- Solve the system of equations and compute reactions ------
[a]=solveq(K,f,bc); %return node displacement
Ed=extract(Edof,a); %write disp in element form
%----- Postprocess ----------------------------------------------
%deformation
num=0;
for i = 1:(m_width + 1)
    for j = 1:(n_length + 1)
        num = num + 1;
        x(i, j) = (j - 1) * uni;
        y(i, j) = (i - 1) * uni;
        z(i, j) = abs(a((num - 1) * 3 + 1)); % Deflection
        theta_x = a((num - 1) * 3 + 2); % Rotation about x-axis
        theta_y = a((num - 1) * 3 + 3); % Rotation about y-axis
        results(num, :) = [num, z(i, j), theta_x, theta_y]; % Store results
        XX(num) = x(i, j);
        YY(num) = y(i, j);     
    end
end

figure
surf(x,y,z)
title('Displacement')

%internel force
for i=1:nie
    F_node(i,:)=Ke(:,:,i)*Ed(i,:)'; % element internel force
end
%plot internal force with patch(matlab commmand)
Shear=zeros((m_width+1)*(n_length+1),1);
Mx=zeros((m_width+1)*(n_length+1),1);
My=zeros((m_width+1)*(n_length+1),1);
right_line_node=[(n_length+1):(n_length+1):m_width*(n_length+1)];
bottom_line_node=[(m_width*(n_length+1)+1):((m_width+1)*(n_length+1)-1)];
last_node=(m_width+1)*(n_length+1);
for i=1:nie
        Shear(Enode(i,2))=F_node(i,1);
        Mx(Enode(i,2))=F_node(i,2);
        My(Enode(i,2))=F_node(i,3);
end
for i=1:length(right_line_node)
    Shear(right_line_node(i))=-F_node(n_length+(i-1)*n_length,4);
    Mx(right_line_node(i))=F_node(n_length+(i-1)*n_length,5);
    My(right_line_node(i))=F_node(n_length+(i-1)*n_length,6);
end
for i=1:length(bottom_line_node)
    Shear(bottom_line_node(i))=-F_node(n_length+(m_width-2)*n_length+i,10);
    Mx(bottom_line_node(i))=-F_node(n_length+(m_width-2)*n_length+i,11);
    My(bottom_line_node(i))=F_node(n_length+(m_width-2)*n_length+i,12);

end
Shear(last_node)=F_node(nie,7);
Mx(last_node)=F_node(nie,8);
My(last_node)=F_node(nie,9);
    
%1.shear plot        
figure
title('Shear')
patch('Faces',Enode(:,2:5),'Vertices',[XX',YY'],'facevertexCdata',Shear,'edgecolor','none','facecolor','interp');%'interp' or 'flat';interp smooth the color
colorbar;
%2.Mx plot
figure
title('Mx')
patch('Faces',Enode(:,2:5),'Vertices',[XX',YY'],'facevertexCdata',Mx,'edgecolor','none','facecolor','interp');%'interp' or 'flat';interp smooth the color
colorbar;
figure
title('My')
patch('Faces',Enode(:,2:5),'Vertices',[XX',YY'],'facevertexCdata',My,'edgecolor','none','facecolor','interp');%'interp' or 'flat';interp smooth the color
colorbar;

%3.scatter
figure;
scatter3(x, y, z,"*");
title('Final Displacement of Each Point');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Displacement (mm)');
colorbar;

% 3.scatter
figure;
scatter3(x, y, zeros(size(z)), '*');
title('Initial Position of Each Point');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Initial Position');
colorbar;

figure;
imagesc(x(1,:), y(:,1), z);
title('Deflection Heatmap');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
colorbar;
axis equal;

figure;
uitable('Data', results, 'ColumnName', {'Point Number', 'Deflection', 'Theta_x', 'Theta_y'}, 'Units', 'Normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
title('Node Deflection and Rotation');

%visualize
figure;
h = surf(x, y, zeros(size(z)));
title('Dynamic Displacement Process');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Displacement (mm)');
colorbar;
zlim([0 max(z(:))]);
% update
num_steps = 100;
for step = 1:num_steps
    current_z = (step / num_steps) * z;
    set(h, 'ZData', current_z);
    drawnow;
    pause(0.1); % gap
end

% visualize result
set(h, 'ZData', z);
drawnow;
% %------------------------ end -----------------------------------








