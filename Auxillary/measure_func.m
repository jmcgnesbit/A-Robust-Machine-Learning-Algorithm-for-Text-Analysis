function measures = measure_func(Theta, meetings)
    
    num_met = size(unique(meetings), 1);
    D = size(meetings, 1);
    measures = zeros(D, 4);
    % For each meeting
    for i = 1:num_met
        index = (meetings == i);
        sub_Theta = Theta(:, index);
        % Computes the mean for each topic across the speakers in meeting i
        pibar = mean(sub_Theta, 2); 
        
        % Bhattacharyya coefficient
        measures(index, 1) = sum(sqrt((sub_Theta + eps) .* (pibar + eps)), 1)';
        % Dot product
        measures(index, 3) = sum(sub_Theta .* pibar, 1)';
        %  Kullback-Leibler
        measures(index, 4) = exp( - sum(pibar .* (log(pibar + eps) - log(sub_Theta + eps)),1));
    end
    
    % Hellinger distance 
    measures(:, 2)  = real(sqrt(1 - measures(:, 1) + 5*eps));
    
end
