clear; clc;

% 定义板的长度和宽度
a = 1; % 板的长度
b = 1; % 板的宽度
E = 2.1e11; % 弹性模量
v = 0.3; % 泊松比
q = 1e6; % 均布荷载

% 定义解析解函数
analytical_solution = @(x, y, D) sum_series(x, y, a, b, q, D);

% 定义厚度值
thickness_values = [0.01];
n_values = linspace(8, 64, 8); % 网格密度
errors = zeros(length(n_values), length(thickness_values));

% 计算不同厚度和网格密度下的误差
for t_idx = 1:length(thickness_values)
    t = thickness_values(t_idx);
    D = E * t^3 / (12 * (1 - v^2)); % 弯曲刚度

    for idx = 1:length(n_values)
        n_length = n_values(idx);
        m_width = n_values(idx);
        uni = 1 / n_length; % 单元大小

        % 初始化有限元模型
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
        ndof = 3;
        Edof = caldof(Enode, ndof); % 自定义函数，假设已定义

        % 生成标签
        num_points = (n_length + 1) * (m_width + 1);
        label = zeros(num_points, 4); % 初始化标签矩阵

        for i = 1:num_points
            label(i, 1) = i; % 点号
            label(i, 2) = 3 * (i - 1) + 1; % 自由度 1
            label(i, 3) = 3 * (i - 1) + 2; % 自由度 2
            label(i, 4) = 3 * (i - 1) + 3; % 自由度 3
        end

        % 元素刚度
        D_mat = hooke(E, v); % 自定义函数，假设已定义
        ep = [t]; % 元素厚度
        nie = size(Enode, 1); % 元素数量
        Ke = cell(nie, 1); % 使用单元数组来存储 Ke 矩阵，避免内存问题
        for i = 1:nie
            Ke{i} = platre(ex(i,:), ey(i,:), ep, D_mat); % 自定义函数，假设已定义
        end

        % 组装全局刚度矩阵，使用稀疏矩阵
        K = sparse(max(max(Edof)), max(max(Edof)));
        for i = 1:nie
            K = assem(Edof(i,:), K, Ke{i}); % 自定义函数，假设已定义
        end

        % 载荷向量和边界条件
        f = sparse(max(max(Edof)), 1);	
        f(1:3:end) = -q * uni * uni; % 假设载荷为均布荷载
        num = 0;

        % 找到边界节点
        boundary_node = [];
        for i = 0:m_width 
            boundary_node = [boundary_node, (i) * (n_length + 1) + 1, (i + 1) * (n_length + 1)];
        end
        boundary_node = unique([boundary_node, 1:n_length + 1, (m_width) * (n_length + 1) + 1:(m_width + 1) * (n_length + 1)]);

        % 将边界条件映射到 Edof
        bc = zeros(length(boundary_node), 2);
        for i = 1:length(boundary_node)
            bc(i,:) = [3 * (boundary_node(i) - 1) + 1, 0]; % 简支边界
        end

        % 求解方程组并计算反应力
        [a] = solveq(K, f, bc); % 自定义函数，假设已定义
        Ed = extract(Edof, a); % 自定义函数，假设已定义

        % 计算挠度和误差
        num = 0;
        z = zeros(m_width + 1, n_length + 1);
        w_analytical = zeros(m_width + 1, n_length + 1);
        for i = 1:(m_width + 1)
            for j = 1:(n_length + 1)
                num = num + 1;
                x = (j - 1) * uni;
                y = (i - 1) * uni;
                z(i, j) = abs(a((num - 1) * 3 + 1)); % 挠度
                w_analytical(i, j) = analytical_solution(x, y, D); % 解析解
            end
        end

        % 计算 L2 误差
        error = (z - w_analytical).^2;
        errors(idx, t_idx) = sqrt(sum(error(:))/num_points); % 计算 L2 误差
    end
end

% 绘制误差与网格密度的关系图
hold on;
for t_idx = 1:length(thickness_values)
    plot(n_values, errors(:, t_idx), '*-', 'DisplayName', sprintf('Thickness = %.2f', thickness_values(t_idx)));
end
xlabel('Grid Parameter');
ylabel('Mean Squared Error');
title('MSE on Different Thickness');
legend show;
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