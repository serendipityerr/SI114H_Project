clear;clc;close all;
%----- Topology -------------------------------------------------
% 定义板的参数
a = 1; % 板的长度
b = 1; % 板的宽度
t = 0.01; % 板的厚度
E = 2.1e11; % 弹性模量
v = 0.3; % 泊松比
D = E * t^3 / (12 * (1 - v^2)); % 弯曲刚度
q = 1e6; % 均布荷载


% 网格密度范围
n_values = 128;
%n_values =linspace(32,128,7); % 最大网格密度调低为 200
zmax = zeros(length(n_values), 1);
for idx = 1:length(n_values)
    n_length = n_values(idx);
    m_width = n_values(idx);
    uni = 1 / n_length; % 单元大小
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

    D=hooke(E,v);%stress=D*strain %% Mechanics of Plate and Shell 
    ep=[t];% thickness of element
    nie=size(Enode,1); %number of element
    Ke = cell(nie, 1); % 使用单元数组来存储 Ke 矩阵，避免内存问题
    for i = 1:nie
        Ke{i} = platre(ex(i,:), ey(i,:), ep, D); % 自定义函数，假设已定义
    end

 % 组装全局刚度矩阵，使用稀疏矩阵
    K = sparse(max(max(Edof)), max(max(Edof)));
    for i = 1:nie
        K = assem(Edof(i,:), K, Ke{i}); % 自定义函数，假设已定义
    end

    f = sparse(max(max(Edof)), 1);	
    f(1:3:end) = -q * uni * uni; % 假设载荷为均布荷载
    num=0;
    %find boundary node
    % 找到边界节点
    boundary_node = [];
    for i = 0:m_width 
        boundary_node = [boundary_node, (i) * (n_length + 1) + 1, (i + 1) * (n_length + 1)];
    end
    boundary_node = unique([boundary_node, 1:n_length + 1, (m_width) * (n_length + 1) + 1:(m_width + 1) * (n_length + 1)]);
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
        end
    end
    zmax(idx) = max(max(z));  %max value of z
end
figure;
plot(n_values, zmax, '-o');
xlabel('q');
ylabel('z');
title('z 与 q 的关系');
grid on;

