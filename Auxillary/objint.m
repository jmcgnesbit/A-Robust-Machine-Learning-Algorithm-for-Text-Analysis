function [beta, beta_nc, measures] = objint(Theta, table_design, index, weights)
    
Theta_speak_meet = [];
for i = 1:size(weights,2)
    Theta_speak_meet = [Theta_speak_meet, (Theta(:, weights{2,i}) * weights{1,i}')];
end




    % Compute measures: [bhat, hell, dotp, klsim]
    measures = measure_func(Theta_speak_meet, table_design.MeetingNumber);
    % Trim to size
    measure_trim = measures(index, :);
    table_design_trim = table_design(index, :);
    
    beta = zeros(4,4);
    beta_nc = zeros(4,4);
    %resid = zeros(sum(index_speaksub), 4);
    
    for i = 1:4
        % Add measure
        table = addvars(table_design_trim, measure_trim(:, i));
        table.Properties.VariableNames(end) = {'Measure'};
        
        % Model with controls
        lmmod = fitlm(table, 'Measure~Transparency+Recession+EPU+Twoday+Stems+Members+Phds');
        beta(i, :) = table2array(lmmod.Coefficients(2,:));
    
        % Model with no controls
        lmmod_nc = fitlm(table, 'Measure~Transparency');
        beta_nc(i,:) = table2array(lmmod_nc.Coefficients(2,:));
            
        % Regress controls without constant or transparency
        lmmod_2 = fitlm(table, 'Measure~Recession+EPU+Twoday+Stems+Members', 'Intercept', false);
        % Save resids
        %resid(:, i) = lmmod_2.Residuals.Raw;
    end
    
end