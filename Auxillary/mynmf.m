function [B, Theta, ii, delta] = mynmf(P, W, K, maxiter, tol)

[V, D] = size(P);

% Initialize B
B = rand(V, K);
B = B ./ sum(B);
% Initialize Theta
Theta = rand(K, D);
Theta = Theta ./ sum(Theta);

% Precompute
WP = W .* P;

for ii = 1: maxiter
    Theta_new = (Theta ./ (B' * W + eps)) .* (B' * ( WP ./ (B * Theta + eps)));
    B_new = (B ./ (W * Theta_new' + eps)) .* (( WP ./ (B * Theta_new + eps))*Theta_new');

    dTheta = max(max(abs(Theta_new - Theta) / (eps + max(max(abs(Theta))))));
    dB = max(max(abs(B_new - B) / (eps+max(max(abs(B))))));
    delta = max(dTheta,dB);
    
    % Check for convergence
    if delta <= tol
        break;
    end
    
    % Update
    Theta = Theta_new;
    B = B_new;  
 
end

% Normalize
B = B ./ sum(B);
Theta = Theta ./ sum(Theta);
end
