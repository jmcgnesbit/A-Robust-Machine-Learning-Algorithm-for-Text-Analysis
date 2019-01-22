function [beta, beta_nc, measures] = objint(Theta, table_design, index)
    

    design_members= [];
design_meeting = [];
Theta_speak_meet = [];
design_stems = [];
for j = 1:num_meeting
    meeting_ind = (dates_cat == j); % find indices of meeting j
    meet_mem = unique(members(meeting_ind)); % find members attending meeting j
    for i = 1:size(meet_mem)
        design_meeting = [design_meeting; j];
        design_members = [design_members; meet_mem(i)];
        member_ind = (members == meet_mem(i));
        combo_ind = and(meeting_ind, member_ind);
        mem_N_d = N_d(combo_ind);
        weights = mem_N_d ./ sum(mem_N_d);
        Theta_speak_meet = [Theta_speak_meet; (weights * store_Theta(1,:, combo_ind)')];
        
        design_stems = [design_stems; sum(any(tdmat(:,combo_ind) ~=0, 2))]; % find non-zero entries in rows (terms) of tdmat
    end
end

design_stems = full(design_stems);




    % Compute measures: [bhat, hell, dotp, klsim]
    measures = measure_func(Theta, table_design.MeetingNumber);
    % Trim to size
    measure_trim = measures(index, :);
    table_design_trim = table_design(index, :);
    
    beta = zeros(4,4);
    beta_nc = zeros(4,4);
    %resid = zeros(sum(index_speaksub), 4);
    
    for i = 1:4
        % Add measure
        table = addvars(table_design_trim, measure_trim(:, i));
        table.Properties.VariableNames(8) = {'Measure'};
        
        % Model with controls
        lmmod = fitlm(table, 'Measure~Transparency+Recession+EPU+Twoday+Stems+Members');
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