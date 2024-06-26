function [K_c, M_c, F_c, nctot] = Constraints(nc, K, M, F)
    % 四边简支边缘条件
    
    % Initialize nctot array to store all constrained DOF
    nctot = zeros(1, numel(nc)); % Each node has 1 constrained DOF for simple support: vertical displacement
    
    % Loop through each constrained node
    for i = 1:numel(nc)
        % For each node, store the DOFs corresponding to the node
        % Each node has 1 DOF for simple support: (i*3-2) for vertical displacement
        nctot(i) = nc(i) * 3 - 2; % constrained DOF
    end
    
    % Remove rows and columns corresponding to constrained DOFs from the stiffness matrix
    K(nctot, :) = [];
    K(:, nctot) = [];
    K_c = K;
    
    % Remove rows and columns corresponding to constrained DOFs from the mass matrix
    M(nctot, :) = [];
    M(:, nctot) = [];
    M_c = M;
    
    % Remove entries corresponding to constrained DOFs from the force vector
    F(nctot) = [];
    F_c = F;
    
    % Return the modified matrices and the list of constrained DOFs
    return
end
