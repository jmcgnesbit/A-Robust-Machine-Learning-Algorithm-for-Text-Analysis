function ydif = rwcyc(y)
% Matlab code to calculate cyclical component based on 2-year-ahead forecast error for baseline case of random walk as recommended in
%      James D. Hamilton, "Why You Should Never Use the Hodrick-Prescott Filter"
%      Review of Economics and Statistics, forthcoming
% input:  y = (T x 1) vector of data, tth element is observation for date t
% output  ydif = (T x 1) vector, tth element is cyclical component for date t

T = size(y,1);
ydif = NaN(T,1);
h = 8;    % default for quarterly data and 2-year horizon
ydif(h+1:T,1) = y(h+1:T,1) - y(1:T-h,1);