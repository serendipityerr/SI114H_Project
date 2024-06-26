clear; clc; close all;

% Some constants
E = 2.1 * 10^11;
gamma = 0.3;
h=0.01;
num = linspace(8, 64, 8);
E_values = [8.16e11,1.155e11,2.1e11];
v_values = [0.42,0.34,0.3];
all_errors = zeros(length(num), length(E_values));

for t_idx = 1:length(E_values)
    E = E_values(t_idx);
    v = v_values(t_idx);
    D = (E * h^3) / (12 * (1 - gamma^2));
    
    for idx = 1:length(num)
        n = num(idx);
        xy = generate_grid_points(n); % 自定义函数，生成网格点
        
        dieta_x = 1 / n;
        N = n + 1;
        
        % Calculate xy without bound
        xy_nobound = zeros((N-2)*(N-2), 2);
        cnt = 0;
        for i = 1:N^2
            col = mod(i, N);
            row = (i - col) / N + 1;
            if (col == 1) || (col == 0) || (row == 1) || (row == N)
                continue;
            else
                cnt = cnt + 1;
                xy_nobound(cnt, 1) = xy(i, 1);
                xy_nobound(cnt, 2) = xy(i, 2);
            end
        end
        
        % Calculate matrix A
        A_row = zeros(1, N-2);
        A_row(1) = -8; A_row(2) = 2;
        A = toeplitz(A_row);
        
        % Calculate matrix I
        I = eye(N-2, N-2);
        
        % Calculate matrix B18
        B18_row = zeros(1, N-2);
        B18_row(1) = 19; B18_row(2) = -8; B18_row(3) = 1;
        B18 = toeplitz(B18_row);
        B18(1, 1) = 18; B18(N-2, N-2) = 18;
        
        % Calculate matrix B19
        B19 = B18 + I;
        
        % Calculate lhs matrix K
        K = [];
        for i = 1:N-2
            K_rows = [];
            if i == 1
                K_rows = [B18 A I zeros(N-2, (N-2)*(N-5))];
            elseif i == 2
                K_rows = [A B19 A I zeros(N-2, (N-2)*(N-6))];
            elseif i == N-2
                K_rows = [zeros(N-2, (N-2)*(N-5)) I A B18];
            elseif i == N-3
                K_rows = [zeros(N-2, (N-2)*(N-6)) I A B19 A];
            elseif (i >= 3) && (i <= N-4)
                K_rows = [zeros(N-2, (N-2)*(i-3)) I A B19 A I zeros(N-2, (N-2)*(N-i-4))];
            end
            K = [K; K_rows];
        end
        
        % Define q(x,y)
        q = (10^6 * dieta_x^4 / D) * ones((N-2)*(N-2), 1);
        
        % Calculate rhs vector f
        f = q;
        
        % Solve equation
        u = K \ f;
        
        % Calculate exact solution
        u_exact = zeros((N-2)*(N-2), 1);
        for k = 1:(N-2)*(N-2)
            for i = 1:50
                for j = 1:50
                    u_exact(k) = u_exact(k) + ...
                        (10^6 * 4 / (pi^6 * D)) * ...
                        (sin(i * pi * xy_nobound(k, 1)) * sin(j * pi * xy_nobound(k, 2)) * ...
                        (1 - (-1)^i) * (1 - (-1)^j)) / ...
                        (i * j * (i^2 + j^2)^2);
                end
            end
        end
        
        % Calculate error
        error = 0.0;
        for k = 1:(N-2)*(N-2)
            error = error + (u(k) - u_exact(k))^2;
        end
        error = sqrt(error / (N-2)^2);
        
        all_errors(idx, t_idx) = error;
    end
end

% 绘制误差与网格密度的关系图
figure;
hold on;
for t_idx = 1:length(E_values)
    plot(num, all_errors(:, t_idx), '*-', 'DisplayName', sprintf('{%.3e,%.2f}',E_values(t_idx), v_values(t_idx)));
end
xlabel('Grid Parameter');
ylabel('Mean Squared Error');
title('MSE on Different E and Gamma');
legend show;
grid on;

