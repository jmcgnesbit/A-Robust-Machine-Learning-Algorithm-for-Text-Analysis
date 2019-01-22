function [upperbound, lowerbound] = optimaldraws(d, epsilon, delta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimaldraws Compute the optimal number of draws from the inside using
    %   Theorem 1.
    %
    % Inputs
    % d - dimension of parameter regin
    % epsilon
    % delta
    %
    % Outputs 
    % upperbound 
    % lowerbound
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lowerbound = (1- epsilon)/(epsilon)* log(1/delta);
    upperbound = ((2*d)/epsilon) * log((2*d)/delta); 

    lowerbound = ceil(lowerbound);
    upperbound = ceil(upperbound);

end