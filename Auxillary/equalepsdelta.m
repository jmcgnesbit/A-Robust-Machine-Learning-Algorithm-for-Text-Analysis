function output = equalepsdelta(dim, M)

output = fzero(@(delta) epsdelta(dim, delta, M), 0.1);
end

function zero = epsdelta(dim, epsilon, M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % epsdelta Function for fzero. Solve for iso-draw curve intersection
    %   with 45 degree line (i.e. fix M and d, solve for epsilon = delta
    %   such that upper bound of theorem 3 holds)
    %
    % Inputs
    % dim - dimension of parameter region
    % delta 
    % M 
    %
    % Outputs
    % zero - want to be equal to zero
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zero = ((2*dim)/epsilon) * log((2*dim)/epsilon) - M;
end