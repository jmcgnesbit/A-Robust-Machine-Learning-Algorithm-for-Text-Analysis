% Things to talk about re:coding
% We use 100 nmfiters, this may not be enough for KL
% If we used MSE instead of KL with would most likely get quicker results -
%   consider changing to the compute hat{P} using KL and then switch to
%   something faster
% I think that the NMF portion could be made quicker in Julia because
%   one can do the multiply in place command, obviously lots of creating
%   temporary arrays eg. (https://github.com/JuliaStats/NMF.jl), could also
%   compute hat{P} and use off the shelf
% Estimating the speaker - meeting - section distributions. In section 3 on
%   their online appendix, they explain how they generate these
%   distributions. The way it is implemented for the NMF is to weight each
%   document by a speaker in a meeting by the number of spoken words and
%   then sum these weighted topics. 


clc 
clear

%% Options
section_ind = 1;
test = 1; % Testing mode. 1 = low iterations, no parallel loop; 0 = high iterations, parellel loop
graph = 1; % Render figures. 1 = render figures; 0 = don't render

% Directory options
addpath(genpath(pwd)) % include subfolders in working directory
% Picture directory
pic_dir = "./Figures";

% Graphing options
fontsize = 12;
fig_height = 2.5;
fig_width = 6;
fig_fontsize = 12;

kmax = 10; % to plot k most often members

% NNMF options
maxiter = 1000;
tol = 0.000001;
K = 40;
outloop = optimaldraws(1, 0.1, 0.1); % 1060 as D = 1
nmfiters = 100;

% 45 degree line for eps and delta
equalepsdelta(1, 120)


%% Load up and initial data cleaning
% Load for_pepe.csv
[txt_dataraw, txt_textraw] = xlsread('for_pepe.csv'); 

docs = txt_textraw(2:end,4);
% Some manual changes (it was coded as TRUE, but matlab didn't like this)
docs(1465) = {'true'};
docs(6014) = {'true'};
docs(15385) = {'true'};
docs(17038) = {'true'};
docs(18560) = {'true'};
docs(18762) = {'true'};

% % Tokenize documents
% cleanDocuments = tokenizedDocument(docs);
% bow = bagOfWords(cleanDocuments); 
% tdmat = bow.Counts';
% 
% % Wordcloud
% figure('Name', 'Wordcloud')
% wordcloud(bow)
% latex_fig(fig_fontsize, fig_width, fig_height);
% print(gcf, '-depsc2', fullfile(pic_dir, 'wordcloud.eps'))

%% Section subsample

docs_sub = docs(txt_dataraw(:,3) == section_ind);
% Section and members
raw_section = txt_dataraw(2:end,3);
section = raw_section(raw_section == section_ind);
members = string(txt_textraw(2:end,2));
members = members(raw_section == section_ind);

% Tokenize documents
cleanDocuments = tokenizedDocument(docs_sub);
bow = bagOfWords(cleanDocuments); 
tdmat = bow.Counts'; %transpose to get term by document
% Check for columns with all zeros
find(all(tdmat == 0))

% Sizes
[V,D] = size(tdmat);

%% NMF
N_d = sum(tdmat,1);
average_words = full(mean(N_d));
% Average number of words per document
fprintf('The average number of words is: %f', average_words);

% Proportions matrix 
prop_matrix = tdmat ./ N_d;

% Weighting matrix V*D matrix. Each column vector j is just N_j repeated V times
W = repmat(N_d, V, 1);

save(['HPCin' num2str(section_ind) '.mat'], 'prop_matrix', 'W')

outloop = 2;
nmfiters = 2;
store_B = zeros(outloop, V, K);
store_Theta = zeros(outloop, K, sub_D);
iterscomp = zeros(outloop,1);
norm = zeros(outloop,1);

for ii = 1:outloop
    ii
    %Run NNMF
    [B, Theta, iterscomp(ii), norm(ii)] = mynmf(prop_matrix, W, K, nmfiters, tol);
    
    store_B(ii, :, :) = B;
    store_Theta(ii, :, :) = Theta;
  
end

load( )



%% Dates and meeting numbers 
dates = xlsread('dates.xlsx'); %Import dates with days
dates_dt = datetime(string(dates), 'InputFormat', 'yyyyMMdd'); % Convert to datetime
dates_cat = grp2idx(categorical(txt_dataraw(txt_dataraw(:,3) == section_ind,1))); % Convert date vector in txt_dataraw to meeting number
num_meeting = size(unique(dates),1);

%% Resampling thing to generate topics for speaker meeting pairs. I think the 
% best way to do this, is to weight each speaker - meeting - injection 
% by its relative contribution

% Number of stems
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
        Theta_speak_meet = [Theta_speak_meet; (weights * store_Theta(:, combo_ind)')];
        
        design_stems = [design_stems; sum(any(tdmat(:,combo_ind) ~=0, 2))]; % find non-zero entries in rows (terms) of tdmat
    end
end

design_stems = full(design_stems);

%% Regression stuff

% Twoday
covariate_dataraw = xlsread('term_document_top2_all_clean.xlsx', 'covariate');
twoday = covariate_dataraw(:,1);


% Get year and month of txt_dates
[txt_y,txt_m] = ymd(dates_dt);

% Transparency dummy
trans_date = datetime('19931116', 'InputFormat', 'yyyyMMdd');
trans_meeting = find(dates == 19931116); % meeting number
trans = [zeros(trans_meeting,1); ones((num_meeting - trans_meeting), 1)]; 
transmetmem = design_members(design_meeting == trans_meeting);

% NBER recession dummy
[reces_raw, reces_text] = xlsread('USREC.xls', 'A1604:B1825');
reces_dates = datetime(reces_text(:,1), 'InputFormat', 'dd-MMM-yy');
[reces_y,reces_m] = ymd(reces_dates);

% EPU index
[EPU_raw, EPU_text] = xlsread('USEPUINDXD.xls', 'A43:B264');
EPU_dates = datetime(EPU_text(:,1), 'InputFormat', 'dd-MMM-yy');
[EPU_y,EPU_m] = ymd(EPU_dates);

% Dates where coincide
coindates = ismember([reces_y, reces_m], [txt_y,txt_m], 'rows');
% Trim to fit
reces_trim = reces_raw(coindates);
EPU_trim = EPU_raw(coindates);

% Create dated table of covariates
table_cov = timetable(dates_dt, trans, reces_trim, EPU_trim, twoday);
table_cov.Properties.DimensionNames(1) = {'Dates'};
table_cov.Properties.VariableNames = {'Transparency', 'Recession', 'EPU', 'Twoday'};
table_cov;

%% Create design mat

design_size = size(design_meeting,1);

design_transparency = zeros(design_size,1);
design_recession = zeros(design_size,1);
design_EPU = zeros(design_size,1);
design_twoday = zeros(design_size,1);

for i = 1:num_meeting
    meeting_ind = (design_meeting == i);
    design_transparency(meeting_ind ) = trans(i);
    design_recession(meeting_ind ) = reces_trim(i);
    design_EPU(meeting_ind) = EPU_trim(i);
    design_twoday(meeting_ind ) = twoday(i);
end


% Create table
table_design = table(design_meeting, design_transparency, design_recession, design_EPU, design_twoday, design_stems, design_members);
table_design.Properties.VariableNames = {'MeetingNumber', 'Transparency', 'Recession', 'EPU', 'Twoday', 'Stems', 'Members'};
table_design;

% Get dates of meetings in window
start_meet = find(dates == 19891114);
end_meet = find(dates == 19971112) - 1;

start_num = find(table_design.MeetingNumber == start_meet, 1, 'first');
end_num = find(table_design.MeetingNumber == end_meet, 1, 'last');

window_index = zeros(design_size, 1);
window_index(start_num:end_num) = 1;
window_index = logical(window_index);

transmem_index = ismember(design_members, transmetmem);

final_index = and(window_index, transmem_index);
% Preallocate
beta = zeros(outloop, 4, 4); 
beta_nc = zeros(outloop, 4, 4);
measure = zeros(outloop, design_size, 4);
iterscomp = zeros(outloop,1);
norm = zeros(outloop,1);

for ii = 1:outloop
     [beta(ii, :, :), beta_nc(ii, :, :), measure(ii,:, :)] = objint(store_Theta, table_design, final_index);  
end


indices = [13016435, 13025769, 13025770, 13025771, 13025772, 13025773];
%indices = [13016435];
chunksize = 20;
outloop2 = size(indices,2)*chunksize;
% % Preallocate
beta = zeros(outloop2, 4, 4); 
beta_nc = zeros(outloop2, 4, 4);
measure = zeros(outloop2, sub_D, 4);
iterscomp = zeros(outloop2,1);
norm = zeros(outloop2,1);
% 
store_B = zeros(outloop2, V, K);
store_Theta = zeros(outloop2, K, sub_D);
   
for i = 1:size(indices,2)
    tempb = load(['B' num2str(indices(i)) '.mat']);
    store_B(1 + chunksize*(i - 1): chunksize*i,:,:) = tempb.store_B;
    temptheta = load(['Theta' num2str(indices(i)) '.mat']);
    store_Theta(1 + chunksize*(i - 1): chunksize*i,:,:) = temptheta.store_Theta;
    tempobs = load(['obsint' num2str(indices(i)) '.mat']);
    beta(1 + chunksize*(i - 1): chunksize*i,:, :) = tempobs.beta;
    beta_nc(1 + chunksize*(i - 1): chunksize*i,:, :) = tempobs.beta_nc;
    measure(1 + chunksize*(i - 1): chunksize*i,:, :) = tempobs.measure;
    clear tempb temptheta tempobs
end

%% Run LDA
mdl = fitlda(sub_tdmat', K, 'solver', 'cgs',  'Verbose', 1,...
    'LogLikelihoodTolerance',0.000001, 'IterationLimit', 1000, ...
    'InitialTopicConcentration', 50, 'WordConcentration', .025*V);

lda_Theta = mdl.DocumentTopicProbabilities';
lda_B = mdl.TopicWordProbabilities';

[lda_beta, lda_beta_nc, lda_measure] = objint(lda_Theta, sub_table_design, window_index);

% 
% %% Product test
% P1 = squeeze(store_B(1,:,:)) * squeeze(store_Theta(1,:,:));
% P2 = squeeze(store_B(2,:,:)) * squeeze(store_Theta(2,:,:));
%% Min and max
measure1 = 2;
measure1_name = 'HD'
measure2 = 3;
measure2_name = 'DP'

min_beta = min(beta(:, :, 1), [], 1);
max_beta = max(beta(:, :, 1), [], 1);

min_beta_nc = min(beta_nc(:, :, 1), [], 1);
max_beta_nc = max(beta_nc(:, :, 1), [], 1);

% Hell
meas1_tab_nc = {lda_beta_nc(measure1,1); lda_beta_nc(measure1,2); ...
            min_beta_nc(measure1); max_beta_nc(measure1)}; 
meas2_tab = {lda_beta(measure1,1); lda_beta(measure1,2); ...
            min_beta(measure1); max_beta(measure1)}; 
% KLSIM
measure2_tab_nc = {lda_beta_nc(measure2,1); lda_beta_nc(measure2,2); ...
            min_beta_nc(measure2); max_beta_nc(measure2)}; 
measure2_tab = {lda_beta(measure2,1); lda_beta(measure2,2); ...
            min_beta(measure2); max_beta(measure2)}; 

out_tab = table(meas1_tab_nc,meas2_tab, measure2_tab_nc,measure2_tab);
out_tab.Properties.VariableNames = {[measure1_name '_nocont'], measure1_name, [measure2_name '_nocont'], measure2_name};
out_tab.Properties.RowNames = {'dtrans', 'se', 'min', 'max'};  

out_tab

% Pvalues
% Hell
lda_beta(measure1,4)
% Hell_nc
lda_beta_nc(measure1,4)
% KLsim
lda_beta(measure2,4)   
% KLsim_nc
lda_beta_nc(measure2,4)

%% Figures
% Average accross members within meeting
avg = zeros(num_meeting, 4);
for j = 1:4  % for each of the 4 measures
    for i = 1:num_meeting - 1 % for each meeting
        index = (dates_cat(index_section) == i); % the documents pertaining to meeting i
        avg(i, j) = mean(lda_measure(index, j)); %take the average
    end
end

figure('Name', 'avg')
meeting_dates = dates_dt(start_meet:end_meet);
% Measure1
plot(meeting_dates, avg(start_meet:end_meet, measure1), 'DisplayName', ...
      measure1_name, 'LineWidth', 1, 'Color', 'blue');
hold on;
% KLSIM
plot(meeting_dates, avg(start_meet:end_meet, measure2), 'DisplayName', ...
      measure2_name, 'LineWidth', 1, 'Color', 'red')
leg = legend;    
set(leg, 'Interpreter', 'latex', 'Fontsize', 14, 'AutoUpdate','off', 'Location', 'northeast')
ylim([0,1])
% Horizontal  
y1=get(gca,'ylim');
plot([trans_date, trans_date], y1, '--', ...
'Color', 'black','LineWidth', 1, 'HandleVisibility','off')

ylabel('Similarity measures', 'fontsize', 14, 'Interpreter', 'latex');
%xlabel('Meeting date', 'fontsize', 14, 'Interpreter', 'latex');
title('Average similarity', 'fontsize',  14, 'Interpreter', 'latex');
% Resize 12 is font size, 10 is length, 3 is width
latex_fig(fig_fontsize, fig_width, fig_height);
tightfig();
print(gcf, '-depsc2', fullfile(pic_dir, 'diffindiff.eps'))
recessionplot;

%% Min / max of averages

avg_measure = zeros(outloop, num_meeting, 4);
for k = 1:outloop % for each iteration
    for j = 1:4 % for each of the 4 measures
        for i =  1:num_meeting - 1 % for each meeting
            index = (dates_cat(index_section) == i); % the documents pertaining to meeting i
            avg_measure(k, i, j) = mean(measure(k, index, j));
        end
    end
end

min_avg = squeeze(min(avg_measure, [], 1));
max_avg = squeeze(max(avg_measure, [], 1));

figure('Name', 'avg2')
% HD
    plot_min = min_avg(start_meet:end_meet, measure1);
    plot_max = max_avg(start_meet:end_meet, measure1);
    % Min
    plot(meeting_dates, plot_min, 'LineWidth', 1, 'Color', 'black');
    hold on;
    % Max
    plot(meeting_dates, plot_max, 'LineWidth', 1, 'Color', 'black');
     plot(meeting_dates, avg(start_meet:end_meet, measure1), 'DisplayName', ...
      'HD', 'LineWidth', 1, 'Color', 'blue');
    % Fill
    h = fill([meeting_dates', fliplr(meeting_dates')], [plot_min', fliplr(plot_max')], 'r');
    set(h,'facealpha',.2)
    
ylim([0,1])
% Horizontal  
y1=get(gca,'ylim');
plot([trans_date, trans_date], y1, '--', ...
'Color', 'black','LineWidth', 1, 'HandleVisibility','off')

ylabel('Hellinger distance', 'fontsize', 14, 'Interpreter', 'latex');
%xlabel('Meeting date', 'fontsize', 14, 'Interpreter', 'latex');
%title('Average similarity', 'fontsize',  14, 'Interpreter', 'latex');
% Resize 12 is font size, 10 is length, 3 is width
latex_fig(fig_fontsize, fig_width, fig_height);
tightfig();
print(gcf, '-depsc2', fullfile(pic_dir, 'hellrange.eps'))
recessionplot;



figure('Name', 'avg3')
% KL
    plot_min = min_avg(start_meet:end_meet, measure2);
    
    plot_max = max_avg(start_meet:end_meet, measure2);
    % Min
    plot(meeting_dates, plot_min, 'LineWidth', 1, 'Color', 'black');
    hold on;
    % Max
    plot(meeting_dates, plot_max, 'LineWidth', 1, 'Color', 'black');
    
    plot(meeting_dates, avg(start_meet:end_meet, measure2), 'DisplayName', ...
      measure2_name, 'LineWidth', 1, 'Color', 'blue');
  
    % Fill
    h = fill([meeting_dates', fliplr(meeting_dates')], [plot_min', fliplr(plot_max')], 'r');
    set(h,'facealpha',.2)
    
ylim([0,1])
% Horizontal  
y1=get(gca,'ylim');
plot([trans_date, trans_date], y1, '--', ...
'Color', 'black','LineWidth', 1, 'HandleVisibility','off')

ylabel('Dot product similarity', 'fontsize', 14, 'Interpreter', 'latex');
%xlabel('Meeting date', 'fontsize', 14, 'Interpreter', 'latex');
%title('Average similarity', 'fontsize',  14, 'Interpreter', 'latex');
% Resize 12 is font size, 10 is length, 3 is width
latex_fig(fig_fontsize, fig_width, fig_height);
tightfig();
print(gcf, '-depsc2', fullfile(pic_dir, 'dprange.eps'))
recessionplot;

%% Isodraw
iso_config = pic_config(fontsize, 3, 4, pic_dir);
isodraw_plot(1, 120, iso_config, 'isodraw_single');
%% Resid

% avg_resid = zeros(num_meeting, 4);
% for j = 1:4
%     for i = 1:num_meeting
%         index = index_meetspeak(index_speaksub, i);
%         avg_resid(i, j) = mean(lda_resid(index, j));
%     end
% end
% 
% figure('Name', 'avg_resid')
% 
% % HD
% plot(meeting_dates, avg_resid(start_meet:end_meet,2), 'DisplayName', ...
%           'HD', 'LineWidth', 1, 'Color', 'blue');
% hold on;
% % KLSIM
% plot(meeting_dates, avg_resid(start_meet:end_meet,4), 'DisplayName', ...
%           'KLSim', 'LineWidth', 1, 'Color', 'red')
% leg = legend;    
% set(leg, 'Interpreter', 'latex', 'Fontsize', 14, 'AutoUpdate','off', 'Location', 'northeast')
% ylim([-0.5,.5])
% % Horizontal  
% y1=get(gca,'ylim');
% plot([trans_date, trans_date], y1, '--', ...
% 'Color', 'black','LineWidth', 1, 'HandleVisibility','off')
% 
% ylabel('Residuals', 'fontsize', 14, 'Interpreter', 'latex');
% %xlabel('Meeting date', 'fontsize', 14, 'Interpreter', 'latex');
% title('Average similarity', 'fontsize',  14, 'Interpreter', 'latex');
% % Resize 12 is font size, 10 is length, 3 is width
% latex_fig(fig_fontsize, fig_width, fig_height);
% tightfig();
% print(gcf, '-depsc2', fullfile(pic_dir, 'resids.eps'))
% recessionplot;


%%   

   
   
% for i = 1:num_met2
%    index = 1:D;
%    index = index(date_count == i);
%    beta_min_avg = mean(bhat_min(index));
%    beta_max_avg = mean(bhat_max(index));
%    beta_lda_avg = mean(lda_bhat(index);
% end
% % Bhat sim


% figure('Name','bhat')
% bhat_min = min(store_bhat, [], 2); 
% bhat_max = max(store_bhat, [], 2);
% 
% % Min
% plot(1:D, bhat_min, 'LineWidth', 1, 'Color', 'black')
% hold on;
% % Max
% plot(1:D, bhat_max, 'LineWidth', 1, 'Color', 'black')
% % Fill
% h = fill([1:D, fliplr(1:D)], [bhat_min', fliplr(bhat_max')], 'r');
% set(h,'facealpha',.2)
% % LDA
% plot(1:D, lda_bhat, 'LineWidth', 1, 'Color', 'blue')
% 
% ylabel('Bhattacharyya similarity', 'fontsize', 12, 'Interpreter', 'latex');
% xlabel('Meeting number', 'fontsize', 12, 'Interpreter', 'latex');
% latex_fig(fontsize, fig_width, fig_height);
% tightfig();
% saveas(gcf, fullfile(pic_dir, 'bhat.png'));
% 
% 
% 
% figure('Name','dotp')
% % Dotp sim
% dotp_min = min(store_dotp, [], 2); 
% dotp_max = max(store_dotp, [], 2);
% 
% % Min
% plot(1:D, dotp_min, 'LineWidth', 1, 'Color', 'black')
% hold on;
% % Max
% plot(1:D, dotp_max, 'LineWidth', 1, 'Color', 'black')
% % Fill
% h = fill([1:D, fliplr(1:D)], [dotp_min', fliplr(dotp_max')], 'r');
% set(h,'facealpha',.2)
% % LDA
% plot(1:D, lda_dotp, 'LineWidth', 1, 'Color', 'black')
% 
% ylabel('Dot product similarity', 'fontsize', 12, 'Interpreter', 'latex');
% xlabel('Meeting number', 'fontsize', 12, 'Interpreter', 'latex');
% latex_fig(fontsize, fig_width, fig_height);
% tightfig();
% 
% 
% 
% % KL sim
% figure('Name','klsim')
% klsim_min = min(store_klsim, [], 2); 
% klsim_max = max(store_klsim, [], 2);
% 
% % Min
% plot(1:D, klsim_min, 'LineWidth', 1, 'Color', 'black')
% hold on;
% % Max
% plot(1:D, klsim_max, 'LineWidth', 1, 'Color', 'black')
% % Fill
% h = fill([1:D, fliplr(1:D)], [klsim_min', fliplr(klsim_max')], 'r');
% set(h,'facealpha',.2)
% % LDA
% plot(1:D, lda_klsim, 'LineWidth', 1, 'Color', 'black')
% 
% ylabel('KL divergence', 'fontsize', 12, 'Interpreter', 'latex');
% xlabel('Meeting number', 'fontsize', 12, 'Interpreter', 'latex');
% latex_fig(fontsize, fig_width, fig_height);
% tightfig();
% saveas(gcf, fullfile(pic_dir, 'klsim.png'));
% 
% 
% 
% %% 
% %% 
% min_t = min(store_Theta, [], 3);
% max_t = max(store_Theta, [], 3);
% min_b = min(store_B, [], 3);
% max_b = max(store_B, [], 3);
% 
% 
% % Plot Theta and Beta across documents for a particular topic
% topic = 1;
% 
% % Theta
% figure('Name','Theta')
% 
% % LDA scatter
% scatter(1:D, lda_Theta(topic,:))
% hold on;
% % Min
% plot(1:D, min_t(topic,:), 'LineWidth', 1, 'Color', 'black')
% hold on;
% % Max
% plot(1:D, max_t(topic,:), 'LineWidth', 1, 'Color', 'black')
% % Fill
% h = fill([1:D, fliplr(1:D)], [min_t(topic,:), fliplr(max_t(topic,:))], 'r');
% set(h,'facealpha',.2)
% 
% ylabel('Theta', 'fontsize', 12, 'Interpreter', 'latex');
% xlabel('Meeting number', 'fontsize', 12, 'Interpreter', 'latex');
% latex_fig(fontsize, fig_width, fig_height);
% tightfig();
% saveas(gcf, fullfile(pic_dir, 'bhat.png'));
% 
% 
% 
% % B
% figure('Name','B')
% % LDA scatter
% scatter(1:V, lda_B(topic,:)')
% hold on;
% % Min
% plot(1:V, min_b(:,topic), 'LineWidth', 1, 'Color', 'black')
% % Max
% plot(1:V, max_b(:,topic), 'LineWidth', 1, 'Color', 'black')
% % Fill
% h = fill([1:V, fliplr(1:V)], [min_b(:,topic)', fliplr(max_b(:,topic)')], 'r');
% set(h,'facealpha',.2)
% 
% ylabel('Theta', 'fontsize', 12, 'Interpreter', 'latex');
% xlabel('Meeting number', 'fontsize', 12, 'Interpreter', 'latex');
% latex_fig(fontsize, fig_width, fig_height);
% tightfig();
% saveas(gcf, fullfile(pic_dir, 'bhat.png'));





%mymembers = ["GREENSPAN"; "KOHN"];
% num_mem = size(mymembers, 1);
% % Check if both members in meeting
% include_meeting = zeros(D,1);
% for i = 1:num_met
%         index = 1:size(txt_data,2);
%         index = index(date_count == i);
%         temp = true;
%         for j = 1:size(mymembers)
%             temp = temp * any(txt_text(index, 2) == mymembers(j)); 
%         end
%         include_meeting(index) = temp;
% end 
% % Chop off irrelevant meetings
% txt_data = txt_data(:, include_meeting == 1);
% txt_text = txt_text((include_meeting == 1), :);
% date_count = date_count(include_meeting == 1);

% meeting_number = unique(meeting_number);
%num_met2 = size(meeting_number,1);
