clear;close all;clc
%---------------------GEOMETRY-----------------------
a=1;
b=1;
q=1e6;
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
D=material.E * material.t^3 / (12 * (1 - material.v^2));

% 定义解析解函数
analytical_solution = @(x, y) sum_series(x, y, a, b, q, D);
%---------------------MESH---------------------------
ex_values = linspace(8,64,15);
ey_values = linspace(8,64,15);
errors = zeros(length(ex_values), 1);
for idx = 1:length(ex_values)
    ex = ex_values(idx);
    ey = ey_values(idx);
    mesh=MeshGenerator(a,b,ex,ey);

    K=StiffnessMatrix(material,mesh,quadrature,ShapeOption);
    M=MassMatrix(material,mesh,quadrature,ShapeOption);
    %---------------------LOADS----------------------------
    
    F=ForceVector(mesh,quadrature,ShapeOption,q);
    % c0=input('\n集中荷载的大小: ');
    % cl=input('\nLato da caricare: ');
    % F=ForceVector(cl,c0,mesh);Lato da caricare
    %---------------------CONSTRAINTS-----------------------
    %nc is nodes of constrained elemment
    nc=[mesh.lato1 mesh.lato2(2:end) mesh.lato3(2:end) mesh.lato4(2:end-1)];%%四边固定
    [K_c,M_c,F_c,nctot]=Constraints(nc,K,M,F);
    %---------------------SOLUTION---------------------------
    [w, ~, ~] = StaticSolver(K_c, F_c, mesh, nctot);

      
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
figure;
plot(ex_values, errors/8, '*-r');
xlabel('Grid parameter');
ylabel('Mean Square Error');
title('MSE on Different Grid Density after Shear Locking');
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
