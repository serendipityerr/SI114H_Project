square_domain
load square_grid.mat
% Some constants
E = 2.1*10^11;
h = 0.01;
gamma = 0.3;
D = (E*h^3)/(12*(1-gamma^2));
dieta_x = 1/n;
N = n+1;
Ns = 10;
d_m = 0.05;

% calculate xy without bound
xy_nobound = zeros((N-2)*(N-2),2);
cnt = 0;
for i = 1:N^2
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

% calculate the nearest Ns points: w_ij, h, k
dist = [];
weight = [];
dist_index = [];
dieta_x = [];
dieta_y = [];
for i = 1:(N-2)^2
    temp_dist = [];
    for j = 1:(N-2)^2
        d = sqrt((xy_nobound(i,1) - xy_nobound(j,1))^2 + (xy_nobound(i,2) - xy_nobound(j,2))^2);
        temp_dist = [temp_dist d];
    end
    [near_temp_dist, near_temp_index] = mink(temp_dist,Ns+1);
    temp_weight = zeros(1,Ns);
    % size(dist) will be 49x6
    dist = [dist; near_temp_dist];
    dist_index = [dist_index; near_temp_index];
    for k = 2:Ns+1
        if near_temp_dist(k) <= d_m
            temp_weight(k-1) = 1 - 6*(near_temp_dist(k)/d_m)^2 + 8*(near_temp_dist(k)/d_m)^3 - 3*(near_temp_dist(k)/d_m)^4;
        end
    end
    weight = [weight; temp_weight];
    temp_dieta_x = [];
    temp_dieta_y = [];
    for m = 2:Ns+1
        temp_dieta_x = [temp_dieta_x (xy_nobound(i,1) - xy_nobound(near_temp_index(m),1))];
        temp_dieta_y = [temp_dieta_y (xy_nobound(i,2) - xy_nobound(near_temp_index(m),2))];
    end
    dieta_x = [dieta_x; temp_dieta_x];
    dieta_y = [dieta_y; temp_dieta_y];
end

% calculate A_P
for i = 1:(N-2)^2
    A_P = zeros(14,14);
    D_P = zeros(14,2);
    coeff = [];
    for j = 1:Ns
        temp_e = zeros(14,1);
        temp_e(1) = dieta_x(i,j);
        temp_e(2) = dieta_y(i,j);
        temp_e(3) = (1/2)*dieta_x(i,j)^2;
        temp_e(4) = dieta_x(i,j)*dieta_y(i,j);
        temp_e(5) = (1/2)*dieta_y(i,j)^2;
        temp_e(6) = (1/6)*dieta_x(i,j)^3;
        temp_e(7) = (1/2)*dieta_x(i,j)^2*dieta_y(i,j);
        temp_e(8) = (1/2)*dieta_x(i,j)*dieta_y(i,j)^2;
        temp_e(9) = (1/6)*dieta_y(i,j)^3;
        temp_e(10) = (1/24)*dieta_x(i,j)^4;
        temp_e(11) = (1/6)*dieta_x(i,j)^3*dieta_y(i,j);
        temp_e(12) = (1/4)*dieta_x(i,j)^2*dieta_y(i,j)^2;
        temp_e(13) = (1/6)*dieta_x(i,j)*dieta_y(i,j)^3;
        temp_e(14) = (1/24)*dieta_y(i,j)^4;
        temp_A = weight(i,j)^2 * (temp_e * temp_e');
        A_P = A_P + temp_A;
        coeff = [coeff, weight(i,j)^2 * temp_e];
    end

    temp_e = zeros(14,1);
    for j = 1:Ns
        temp_e(1) = temp_e(1) + weight(i,j)^2 * dieta_x(i,j);
        temp_e(2) = temp_e(2) + weight(i,j)^2 * dieta_y(i,j);
        temp_e(3) = temp_e(3) + weight(i,j)^2 * (1/2)*dieta_x(i,j)^2;
        temp_e(4) = temp_e(4) + weight(i,j)^2 * dieta_x(i,j)*dieta_y(i,j);
        temp_e(5) = temp_e(5) + weight(i,j)^2 * (1/2)*dieta_y(i,j)^2;
        temp_e(6) = temp_e(6) + weight(i,j)^2 * (1/6)*dieta_x(i,j)^3;
        temp_e(7) = temp_e(7) + weight(i,j)^2 * (1/2)*dieta_x(i,j)^2*dieta_y(i,j);
        temp_e(8) = temp_e(8) + weight(i,j)^2 * (1/2)*dieta_x(i,j)*dieta_y(i,j)^2;
        temp_e(9) = temp_e(9) + weight(i,j)^2 * (1/6)*dieta_y(i,j)^3;
        temp_e(10) = temp_e(10) + weight(i,j)^2 * (1/24)*dieta_x(i,j)^4;
        temp_e(11) = temp_e(11) + weight(i,j)^2 * (1/6)*dieta_x(i,j)^3*dieta_y(i,j);
        temp_e(12) = temp_e(12) + weight(i,j)^2 * (1/4)*dieta_x(i,j)^2*dieta_y(i,j)^2;
        temp_e(13) = temp_e(13) + weight(i,j)^2 * (1/6)*dieta_x(i,j)*dieta_y(i,j)^3;
        temp_e(14) = temp_e(14) + weight(i,j)^2 * (1/24)*dieta_y(i,j)^4;
    end
    
    u = [];
    for j = 0:Ns
        if j == 0
            sol = A_P \ temp_e;
            u = [u,(-1)*sol];
        end
        if j >= 1
            sol = A_P \ coeff(:,j);
            u = [u,sol];
        end
    end

end
