clc 
clear

%% Options

test = 0; % Testing mode. 1 = low iterations, no parallel loop; 0 = high iterations, parellel loop
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
K = 6;
if test == 1
    outloop = 10;
else
    outloop = optimaldraws(1, 0.01, .01); % 1060 as D = 1
end

% Members to plot and regress
mymembers = ["GREENSPAN"; "STERN"; "KOHN"; "PARRY"];

%% Load up and initial data cleaning
% Load up the term document matrix from Barry
[txt_dataraw, txt_textraw] = xlsread('term_document_top2_all_clean.xlsx', 'Sheet1'); 
txt_dataraw = txt_dataraw(2:end, 2:end)'; % Transpose to get into V * D and clip off headings 
txt_textraw = string(txt_textraw(2:end, :)); % Trim heading
% Copy to working variable
txt_data = txt_dataraw;
txt_text = txt_textraw;

% Drop the columns of zeros (no words for speaker)
zerocol = find(all(txt_data==0,1));
txt_data(:, zerocol) = [];
% Drop for text as well
txt_text(zerocol, :) = [];

% Sizes
[V,D] = size(txt_data);

% Dates and meeting numbers 
date_index = find(~cellfun(@isempty,txt_text(:,1))); % find index of nonempty entries
txt_dates = datetime(txt_text(date_index, 1), 'InputFormat', 'yyyyMMdd'); % Convert to datetime
num_meeting = size(date_index, 1);
date_index = [date_index; D + 1];

% Factor for meeting number
meeting_number = zeros(D,1);
for i = 1:num_meeting 
        meeting_number(date_index(i) : date_index(i + 1) - 1) = i;
end
index_meeting = logical(dummyvar(meeting_number));

% Transparency meeting
trans_meeting = meeting_number(txt_text(:,1) == '19931116'); % meeting number
trans_date = datetime('19931116', 'InputFormat', 'yyyyMMdd'); % meeting in datetime

% Find indexs of the start and end dates
start_d = find(txt_text(:,1) == '19891114');
end_d = find(txt_text(:,1) == '19971112') - 1;


%% Members
members = string(txt_text(:,2));

% There is a member D who should be "LINDSEY"
members(members == "D") = "LINDSEY";

% List of members and dummies
listofmembers = unique(members);
num_members = size(listofmembers, 1);
membersdummies = logical(dummyvar(categorical(members)));

% Take attendence
attendence = zeros(num_meeting, num_members);
for i = 1 : num_meeting
    for j = 1 : num_members
        if any(members(index_meeting(:, i)) == listofmembers(j))
            attendence(i, j) = 1;
        end
    end
end

% Average number of attendees
fprintf('The average number of attendees is: %d', mean(sum(attendence, 2)));

% Find kmax most often members
[~, koftenmeet] = maxk(sum(attendence, 1), kmax);
kmaxlist = listofmembers(koftenmeet);
kmaxattendence = attendence(:, koftenmeet);

if graph == 1
    % Plot members
    figure('Name', 'attendence') 
    for i = 1:kmax
       scatter(txt_dates, i * kmaxattendence(:, i), 'filled')
        hold on;
    end
    % Change ylim
    ylim([.5, kmax + 0.5])
    % Put labels
    yticklabels(kmaxlist)
    % Horizontal  
    y1 = get(gca,'ylim');
    plot([trans_date, trans_date], y1, '--', ...
    'Color', 'black','LineWidth', 1, 'HandleVisibility','off')
    %ylabel('Residuals', 'fontsize', 14, 'Interpreter', 'latex');
    %xlabel('Meeting date', 'fontsize', 14, 'Interpreter', 'latex');
    %title('Residuals', 'fontsize',  14, 'Interpreter', 'latex');
    % Resize 12 is font size, 10 is length, 3 is width
    latex_fig(fig_fontsize, fig_width, fig_height);
    tightfig();
    print(gcf, '-depsc2', fullfile(pic_dir, 'attendence.eps'))
    recessionplot;
end


% Member subsample
% Index 
index_speaker = ismember(members, mymembers);
% Index to subsample
start_d = find(txt_text(:,1) == '19891114');
end_d = find(txt_text(:,1) == '19971112') - 1;
index_subsample = false(D, 1);
index_subsample(start_d:end_d) = 1;

% Joint index for speakers and in subsample 583 index with 252 trues
index_speaksub = index_speaker .* index_subsample;
index_speaksub = logical(index_speaksub(index_speaker));
% Joint index for meeting and subsample
index_meetspeak = index_meeting(index_speaker,:);

% List of included meetings
start_meet = find(date_index == start_d);
end_meet = find(date_index == end_d +1) - 1;

% Average words
mean(sum(txt_data(:, index_speaker)))


%% Regression data
% Get year and month of txt_dates
[txt_y,txt_m] = ymd(txt_dates);

% Transparency dummy
trans = [zeros(trans_meeting,1); ones((num_meeting - trans_meeting), 1)]; 

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

% 2 day meeting (column 1) and stems (column 2) 
covariate_dataraw = xlsread('term_document_top2_all_clean.xlsx', 'covariate');
twoday = covariate_dataraw(:,1);

% Number of stems
stem_raw = zeros(D, 1);
for i = 1:D
    stem_raw(i) = nnz(txt_data(:,i)); % number of nonzero in document 
end

% Create dated table of covariates
table_cov = timetable(txt_dates, trans, reces_trim, EPU_trim, twoday);
table_cov.Properties.DimensionNames(1) = {'Dates'};
table_cov.Properties.VariableNames = {'Transparency', 'Recession', 'EPU', 'Twoday'};

%% Create design mat

transparency = zeros(D ,1);
recession = zeros(D ,1);
EPU = zeros(D ,1);
twoday_des = zeros(D , 1);

for i = 1:num_meeting
    transparency(index_meeting(:,i)) = trans(i);
    recession(index_meeting(:,i)) = reces_trim(i);
    EPU(index_meeting(:,i)) = EPU_trim(i);
    twoday_des(index_meeting(:,i)) = twoday(i);
end

% Create table
table_design = table(meeting_number, transparency, recession, EPU, twoday_des, stem_raw, members);
table_design.Properties.VariableNames = {'MeetingNumber', 'Transparency', 'Recession', 'EPU', 'Twoday', 'Stems', 'Members'};
table_design = table_design(index_speaker, :);

table_design_trim = table_design(index_speaksub, :)
%% NNMF
% Term doc matrix with only relevant speakers
nnmf_data = txt_data(:, index_speaker);
N_d = sum(nnmf_data,1);
% Average number of words per document
fprintf('The average number of words is: %f', mean(N_d));

% Proportions matrix 
prop_matrix = nnmf_data ./ N_d;

% Weighting matrix V*D matrix. Each column vector j is just N_j repeated V times
W = repmat(N_d, V, 1);

% Preallocate
beta = zeros(outloop, 4, 4);
beta_nc = zeros(outloop, 4, 4);
measure = zeros(outloop, sum(index_speaker), 4);

% store_bhat = zeros(D, outloop);
% store_hell = zeros(D, outloop);
% store_dotp = zeros(D, outloop);
% store_klsim = zeros(D, outloop);

for ii = 1:outloop
    ii
    
    % Run NNMF
    [B, Theta] = mynmf(prop_matrix, W, K, 2, maxiter, tol, tol);
    
%     % Turn off storage of B and Theta for the moment
%     store_B(:,:,ii) = B;
%     store_Theta(:,:,ii) = Theta;
    
    [beta(ii, :, :), beta_nc(ii, :, :), measure(ii,:, :), ~] = objint(Theta, table_design, index_speaksub, index_meetspeak);    
end

%% Run LDA
mdl = fitlda(txt_data(:, index_speaker)', K, 'solver', 'cgs',  'Verbose', 1,...
    'LogLikelihoodTolerance',0.00001, 'IterationLimit', 1000, ...
    'InitialTopicConcentration', 50, 'WordConcentration', .025*V);
    
lda_Theta = mdl.DocumentTopicProbabilities';
lda_B = mdl.TopicWordProbabilities';

[lda_beta, lda_beta_nc, lda_measure, lda_resid] = objint(lda_Theta, table_design, index_speaksub, index_meetspeak); 

%% Min and max

min_beta = min(beta(:, :, 1), [], 1);
max_beta = max(beta(:, :, 1), [], 1);

min_beta_nc = min(beta_nc(:, :, 1), [], 1);
max_beta_nc = max(beta_nc(:, :, 1), [], 1);

% Hell
hell_tab_nc = {lda_beta_nc(2,1); lda_beta_nc(2,2); ...
            min_beta_nc(2); max_beta_nc(2)}; 
hell_tab = {lda_beta(2,1); lda_beta(2,2); ...
            min_beta(2); max_beta(2)}; 
% KLSIM
klsim_tab_nc = {lda_beta_nc(4,1); lda_beta_nc(4,2); ...
            min_beta_nc(4); max_beta_nc(4)}; 
klsim_tab = {lda_beta(4,1); lda_beta(4,2); ...
            min_beta(4); max_beta(4)}; 

out_tab = table(hell_tab_nc,hell_tab, klsim_tab_nc,klsim_tab);
out_tab.Properties.VariableNames = {'HD_nocont', 'HD', 'KLSIM_nc', 'KLSIM'};
out_tab.Properties.RowNames = {'dtrans', 'se', 'min', 'max'};  

out_tab

% Pvalues
% Hell
lda_beta(2,4)
% Hell_nc
lda_beta_nc(2,4)
% KLsim
lda_beta(4,4)
% KLsim_nc
lda_beta_nc(4,4)

%% Figures
% Average accross members within meeting
avg = zeros(num_meeting, 4);
for j = 1:4
    for i = 1:num_meeting
        index = index_meetspeak(:, i);
        avg(i, j) = mean(lda_measure(index, j));
    end
end

figure('Name', 'avg')
meeting_dates = txt_dates(start_meet:end_meet);
% HD
plot(meeting_dates, avg(start_meet:end_meet, 2), 'DisplayName', ...
      'HD', 'LineWidth', 1, 'Color', 'blue');
hold on;
% KLSIM
plot(meeting_dates, avg(start_meet:end_meet, 4), 'DisplayName', ...
      'KLSim', 'LineWidth', 1, 'Color', 'red')
leg = legend;    
set(leg, 'Interpreter', 'latex', 'Fontsize', 14, 'AutoUpdate','off', 'Location', 'east')
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
for k = 1:outloop
    for j = 1:4
        for i = 1:num_meeting
            index = index_meetspeak(:, i);
            avg_measure(k, i, j) = mean(measure(k, index, j));
        end
    end
end

min_avg = squeeze(min(avg_measure, [], 1));
max_avg = squeeze(max(avg_measure, [], 1));

figure('Name', 'avg2')
% HD
    plot_min = min_avg(start_meet:end_meet, 2);
    plot_max = max_avg(start_meet:end_meet, 2);
    % Min
    plot(meeting_dates, plot_min, 'LineWidth', 1, 'Color', 'black');
    hold on;
    % Max
    plot(meeting_dates, plot_max, 'LineWidth', 1, 'Color', 'black');
     plot(meeting_dates, avg(start_meet:end_meet, 2), 'DisplayName', ...
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
    plot_min = min_avg(start_meet:end_meet, 4);
    
    plot_max = max_avg(start_meet:end_meet, 4);
    % Min
    plot(meeting_dates, plot_min, 'LineWidth', 1, 'Color', 'black');
    hold on;
    % Max
    plot(meeting_dates, plot_max, 'LineWidth', 1, 'Color', 'black');
    
    plot(meeting_dates, avg(start_meet:end_meet, 4), 'DisplayName', ...
      'HD', 'LineWidth', 1, 'Color', 'blue');
  
    % Fill
    h = fill([meeting_dates', fliplr(meeting_dates')], [plot_min', fliplr(plot_max')], 'r');
    set(h,'facealpha',.2)
    
ylim([0,1])
% Horizontal  
y1=get(gca,'ylim');
plot([trans_date, trans_date], y1, '--', ...
'Color', 'black','LineWidth', 1, 'HandleVisibility','off')

ylabel('KL similarity', 'fontsize', 14, 'Interpreter', 'latex');
%xlabel('Meeting date', 'fontsize', 14, 'Interpreter', 'latex');
%title('Average similarity', 'fontsize',  14, 'Interpreter', 'latex');
% Resize 12 is font size, 10 is length, 3 is width
latex_fig(fig_fontsize, fig_width, fig_height);
tightfig();
print(gcf, '-depsc2', fullfile(pic_dir, 'klrange.eps'))
recessionplot;
%% Resid

avg_resid = zeros(num_meeting, 4);
for j = 1:4
    for i = 1:num_meeting
        index = index_meetspeak(index_speaksub, i);
        avg_resid(i, j) = mean(lda_resid(index, j));
    end
end

figure('Name', 'avg_resid')

% HD
plot(meeting_dates, avg_resid(start_meet:end_meet,2), 'DisplayName', ...
          'HD', 'LineWidth', 1, 'Color', 'blue');
hold on;
% KLSIM
plot(meeting_dates, avg_resid(start_meet:end_meet,4), 'DisplayName', ...
          'KLSim', 'LineWidth', 1, 'Color', 'red')
leg = legend;    
set(leg, 'Interpreter', 'latex', 'Fontsize', 14, 'AutoUpdate','off', 'Location', 'northeast')
ylim([-0.5,.5])
% Horizontal  
y1=get(gca,'ylim');
plot([trans_date, trans_date], y1, '--', ...
'Color', 'black','LineWidth', 1, 'HandleVisibility','off')

ylabel('Residuals', 'fontsize', 14, 'Interpreter', 'latex');
%xlabel('Meeting date', 'fontsize', 14, 'Interpreter', 'latex');
title('Average similarity', 'fontsize',  14, 'Interpreter', 'latex');
% Resize 12 is font size, 10 is length, 3 is width
latex_fig(fig_fontsize, fig_width, fig_height);
tightfig();
print(gcf, '-depsc2', fullfile(pic_dir, 'resids.eps'))
recessionplot;


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
