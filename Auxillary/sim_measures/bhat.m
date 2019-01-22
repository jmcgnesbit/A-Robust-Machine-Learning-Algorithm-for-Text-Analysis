function out = bhat(Theta, date_count)
    
    D = size(Theta,2);
    out = zeros(D, 1);
    num_met = max(date_count);
   
    for i = 1:num_met 
        index = 1:D;
        index = index(date_count == i);
        % Take the mean of the topic composition for each meeting (i.e. for
        % the first 25 documents, take the mean accross the rows).
        pibar = mean(Theta(:, index), 2);
        % Compute the Bhat coefficient in vector form. Multiply each
        % pi_{it} by pibar to create a matrix, take the sqrt, then sum down the rows
        out(index) = sum(sqrt(Theta(:, index) .* pibar), 1)';
    end
    
end
