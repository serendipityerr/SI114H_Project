function xy = generate_grid_points(n)
    %GENERATE_GRID_POINTS Generates n x n grid points in (0,1) x (0,1) coordinate space
    % grid_points = generate_grid_points(n);
    %
    % Input:
    % n: A positive integer representing the number of divisions in each dimension
    %
    % Output:
    % grid_points: A matrix of size n^2 x 2 containing the grid points
    
    if ~isscalar(n) || n <= 0 || mod(n,1) ~= 0
        error('Input must be a positive integer.');
    end
    
    % Generate grid points
    x_vals = linspace(0, 1, n+1);
    y_vals = linspace(0, 1, n+1);
    nvtx=(n+1)*(n+1);
    [X, Y] = meshgrid(x_vals, y_vals);
    xx=reshape(X',nvtx,1);
    yy=reshape(Y',nvtx,1);
    xy=[xx(:),yy(:)];
    
end

