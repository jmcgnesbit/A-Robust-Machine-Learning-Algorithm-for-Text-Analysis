% Matlab code to illustrate calculation of cyclical component based on 2-year ahead forecast error as recommended in 
% Matlab code to calculate cyclical component based on 2-year-ahead forecast error from linear regression as recommended in
%      James D. Hamilton, "Why You Should Never Use the Hodrick-Prescott Filter"
%      Review of Economics and Statistics, forthcoming

% read in data
load employment.csv
    % col 1 = date
    % col 2 = seasonally adjusted nonfarm payrolls
    % col 3 = seasonally unadjusted nonfarm payrolls


yreg = regcyc(100*log(employment(:,2)));
ydif = rwcyc(100*log(employment(:,2)));
disp('Output')
disp('First column is date')
disp('Second column is raw data')
disp('Third column is cyclical component for date t (calculated by regression)')
disp('Fourth column is cyclical component for date t (calculated under assumption of random walk)')
format bank
disp([employment(:,1) 100*log(employment(:,2)) yreg ydif])