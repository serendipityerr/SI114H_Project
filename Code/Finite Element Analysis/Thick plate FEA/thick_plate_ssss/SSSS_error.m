clear all;close all;clc
%---------------------GEOMETRY-----------------------
a=1;
b=1;
q=1e6;
%---------------------ELEMENT------------------------
ShapeOption='Q4';%Q8 Q9
quadrature=GaussQuadrature('gauss1');% Q4 element points and weights
%---------------------MATERIAL------------------------
%SI mm t MPa N
material=struct();
material.t=0.01;%thickness
material.E=2.1*10^11;%2.1e5MPa
material.v=0.3;%Poison ratio
material.rho=2.7e-9;
material.G=material.E/(2+2*material.v);
D=material.E * material.t^3 / (12 * (1 - material.v^2));

% 定义解析解函数
analytical_solution = @(x, y) sum_series(x, y, a, b, q, D);
%---------------------MESH---------------------------
ex_values = linspace(10,60,11);
ey_values = linspace(10,60,11);
errors = zeros(length(ex_values), 1);
for idx = 1:length(ex_values)
    ex = ex_values(idx);
    ey = ey_values(idx);
    mesh=MeshGenerator(a,b,ex,ey);

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
    F=ForceVector(mesh,quadrature,ShapeOption,q);%F
    % F=ForceVector(cl,c0,mesh)
    %---------------------1e-3-----------------------
    %nc is nodes of constrained elemment
    nc=[mesh.lato1 mesh.lato2(2:end) mesh.lato3(2:end) mesh.lato4(2:end-1)];
    [K_c,M_c,F_c,nctot]=Constraints(nc,K,M,F);
    %---------------------SOLUTION---------------------------
    [w,thetax,thethy]=StaticSolver(K_c,F_c,mesh,nctot);
     % 计算 L2 误差
    num = 0;
    z = zeros(ey + 1, ex + 1);
    w_analytical = zeros(ey + 1, ex + 1);
    uni_x = a / ex;
    uni_y = b / ey;
    num_points = (ex+1)^2;
    for i = 1:(ey + 1)
        for j = 1:(ex + 1)
            num = num + 1;
            x = (j - 1) * uni_x;
            y = (i - 1) * uni_y;
            z(i, j) = abs(w(num)); % 挠度
            w_analytical(i, j) = analytical_solution(x, y); % 解析解
        end
    end    
    error = (z - w_analytical).^2;
    errors(idx) = sqrt(sum(error(:))/num_points); % 计算 L2 误差
end
%---------------------PLOT-------------------------------
figure;
plot(ex_values, errors, '*-r');
xlabel('网格密度');
ylabel('L2 误差');
title('不同网格密度下的 L2 误差');
grid on;


% 定义求和函数
function w = sum_series(x, y, a, b, q, D)
    w = 0;
    for m = 1:50
        for n = 1:50
            w = w + (4 * q / (pi^6 * D)) * (sin(m * pi * x / a) * sin(n * pi * y / b)) * (1-(-1)^m) * (1-(-1)^n) / (m * n * ((m^2 / a^2 + n^2 / b^2)^2));
        end
    end
end

