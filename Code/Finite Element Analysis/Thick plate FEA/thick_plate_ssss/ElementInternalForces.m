function [shear, mx, my] = ElementInternalForces(material, mesh, element_nodes, element_displacement)
    % Extract material properties
    E = material.E;
    v = material.v;
    t = material.t;

    % Initialize output forces
    shear = zeros(length(element_nodes), 1);
    mx = zeros(length(element_nodes), 1);
    my = zeros(length(element_nodes), 1);

    % Define shape functions and their derivatives for Q4 element
    N = @(xi, eta) 0.25 * [(1 - xi) * (1 - eta); (1 + xi) * (1 - eta); (1 + xi) * (1 + eta); (1 - xi) * (1 + eta)];
    dN_dxi = @(xi, eta) 0.25 * [-(1 - eta), -(1 - xi); (1 - eta), -(1 + xi); (1 + eta), (1 + xi); -(1 + eta), (1 - xi)];

    % Gauss quadrature points for Q4 element
    gauss_points = [-1, 1, 1, -1; -1, -1, 1, 1] / sqrt(3);

    for i = 1:size(gauss_points, 2)
        xi = gauss_points(1, i);
        eta = gauss_points(2, i);

        % Shape function derivatives at Gauss point
        dN = dN_dxi(xi, eta);

        % Jacobian matrix
        J = dN' * mesh.cordinates(element_nodes, 1:2);
        detJ = det(J);
        invJ = inv(J);

        % Strain-displacement matrix B
        B = zeros(3, 12);
        for j = 1:4
            dNxy = invJ * dN(j, :)';
            B(1, 3*j-2) = dNxy(1); % For w
            B(2, 3*j-1) = dNxy(2); % For theta_x
            B(3, 3*j) = dNxy(1);   % For theta_y
            % Adjust B matrix as per your specific element formulation
        end

        % Ensure element_displacement is a 12x1 vector
        if length(element_displacement) ~= 12
            error('element_displacement must be a 12x1 vector for a Q4 element with 3 DOFs per node.');
        end

        % Calculate strain
        strain = B * element_displacement;

        % Plane stress constitutive matrix
        D = (E / (1 - v^2)) * [1, v, 0; v, 1, 0; 0, 0, (1 - v) / 2];

        % Stress
        stress = D * strain;

        % Calculate forces
        shear(i) = -stress(3) * t; % Shear stress * thickness
        mx(i) = -stress(1) * t;    % Moment in x direction
        my(i) = -stress(2) * t;    % Moment in y direction
    end
end
