% square_domain
% load square_grid.mat

num = [128];
all_error = zeros(length(num),1);
for idx = 1: length(num)
    n = num(idx);
    xy = generate_grid_points(n);
    % Some constants
    E = 2.1*10^11;
    h = 0.01;
    gamma = 0.3;
    D = (E*h^3)/(12*(1-gamma^2));
    dieta_x = 1/n;
    N = n+1;
    
    % calculate xy without bound
    xy_nobound = zeros((N-2)*(N-2),2);
    cnt = 0;
    for i=1:N^2
        col = mod(i,N);
        row = (i - col) / N + 1;
        if (col == 1) || (col == 0) || (row == 1) || (row == N)
            continue;
        else
            cnt = cnt + 1;
            xy_nobound(cnt,1) = xy(i,1);
            xy_nobound(cnt,2) = xy(i,2);
        end
    end
    
    
    % calculate matrix A
    A_row = zeros(1,N-2);
    A_row(1) = -8; A_row(2) = 2;
    A = toeplitz(A_row);
    % calculate matrix I
    I = eye(N-2,N-2); % 单位矩阵
    % calculate matrix B21
    B21_row = zeros(1,N-2);
    B21_row(1) = 21; B21_row(2) = -8; B21_row(3) = 1;
    B21 = toeplitz(B21_row);
    B21(1,1) = 22; B21(N-2,N-2) = 22;
    % calculate matrix B20
    B20 = B21 - I;
    % calculate lhs matrix K
    K = [];
    for i=1:N-2
        K_rows = [];
        if i == 1
            K_rows = [B21 A I zeros(N-2,(N-2)*(N-5))];
        end
        if i == 2
            K_rows = [A B20 A I zeros(N-2,(N-2)*(N-6))];
        end
        if i == N-2
            K_rows = [zeros(N-2,(N-2)*(N-5)) I A B21];
        end
        if i == N-3
            K_rows = [zeros(N-2,(N-2)*(N-6)) I A B20 A];
        end
        if (i >= 3) && (i <= N-4)
            K_rows = [zeros(N-2,(N-2)*(i-3)) I A B20 A I zeros(N-2,(N-2)*(N-i-4))];
        end
        K = [K;K_rows];
    end
    %size(K)
    
    % Define q(x,y)
    q = (10^6 * dieta_x^4 / D)*ones((N-2)*(N-2),1);
    % calculate rhs vector f
    f = q;
    
    % Solve equation
    u = K \ f;
    
    % % calculate exact solution
    % u_exact = zeros((N-2)*(N-2),1);
    % for k = 1 : (N-2)*(N-2)
    %     for i = 1 : 50
    %         for j = 1 : 50
    %             u_exact(k) = u_exact(k) + (sin(i*pi*xy_nobound(k,1))*sin(j*pi*xy_nobound(k,2))*(1-(-1)^i)*(1-(-1)^j))/(i*j*(i^2+j^2)^2);
    %         end
    %     end
    % end
    % 
    % % calculate error
    % error = 0.0;
    % for k = 1 : (N-2)*(N-2)
    %     error = error + (u(k)-u_exact(k))^2;
    % end
    % error = sqrt(error / (N-2)^2)
    
%     u_with_bound = zeros(N^2,1);
%     cnt1 = 0;
%     for i=1:N^2
%         col = mod(i,N);
%         row = (i - col) / N + 1;
%         if (col == 1) || (col == 0) || (row == 1) || (row == N)
%             u_with_bound(i) = 0;
%             continue;
%         else
%             cnt1 = cnt1 + 1;
%             u_with_bound(i) = u(cnt1);
%         end
%     end
%     
%     % plot
%     x_plot = reshape(xy(:,1),N,N);
%     y_plot = reshape(xy(:,2),N,N);
%     u_plot = reshape(u_with_bound,N,N);
%     surf(x_plot,y_plot,u_plot)
end
all_errors



