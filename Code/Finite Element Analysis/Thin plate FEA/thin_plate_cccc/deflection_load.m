clear all;clc;close all;
%----- Topology -------------------------------------------------
%square plate with n element each side
n_length=64;%number of element in long side
m_width=64;%number of element in wide side
uni=1/64;%element size

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
%qs = 10 .^ (1:10);
qs = 10^6 *(1:10);
zmax = zeros(1, 10);
for index = 1:length(qs)  %from 10^1 to 10^10
    q = qs(index);
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
            % theta_x = a((num - 1) * 3 + 2); % Rotation about x-axis
            % theta_y = a((num - 1) * 3 + 3); % Rotation about y-axis
            % results(num, :) = [num, z(i, j), theta_x, theta_y]; % Store results
            % XX(num) = x(i, j);
            % YY(num) = y(i, j);     
        end
    end
    zmax(index) = max(max(z));  %max value of z
end
figure;
plot(qs, zmax, '-o');
xlabel('q');
ylabel('z');
title('z 与 q 的关系');
grid on;

