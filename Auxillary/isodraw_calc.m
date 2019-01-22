function epsilon = isodraw_calc(dim, delta, M)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % isodraw_calc Function to compute iso-draw curves. Upper bound in
    %   theorem 3 rearranged for epsilon
    %
    % Inputs
    % dim - dimension of parameter region
    % delta 
    % M 
    %
    % Outputs 
    % epsilon
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    epsilon = ((2*dim)/M) * log((2*dim)/delta);
end