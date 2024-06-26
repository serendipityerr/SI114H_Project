n_values = linspace(32,128,7);
w_max_software = [0.0642786, 0.0652436, 0.0654653, 0.0657163, 0.0657599, 0.0657677, 0.0657762];
errors_fdm = [1.4598e-04,4.1057e-05,1.8785e-05,1.0703e-05,6.9011e-06,4.8155e-06,3.5498e-06,2.7245e-06,2.1568e-06,1.7497e-06,1.4478e-06,1.2178e-06,1.0385e-06,8.9611e-07,7.8110e-07,6.8689e-07];
errors_fea = [8.9602e-04,2.3510e-04,1.0637e-04,6.0395e-05,3.8875e-05,2.7102e-05,1.9968e-05,1.5320e-05,1.2125e-05,9.8343e-06,8.1364e-06,6.8431e-06,5.8353e-06,5.0349e-06,4.3883e-06,3.8589e-06];
% 绘制误差与网格密度的关系图
figure;
plot(n_values, errors_fdm, '*-b');
xlabel('Grid Parameter');
ylabel('Mean Squared Error');
title('MSE on Different Grid Density');
hold on;
plot(n_values, errors_fea, '*-r');
legend("FDM", "FEA")
grid on;