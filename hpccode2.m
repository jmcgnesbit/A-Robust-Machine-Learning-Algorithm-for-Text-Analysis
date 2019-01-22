function hpccode2(j)

%parpool('local', 'SpmdEnabled', false)

addpath(genpath(pwd)) % include subfolders in working directory

% NNMF options
nmfiters = 100;
tol = 0.01;
K = 40;
outloop = 50;%optimaldraws(1, 0.01, .01); % 1060 as D = 1
chunks = 1;

% Load in data
S = load('HPCin2.mat');

W = S.W;
prop_matrix = S.prop_matrix;
clear S

[V, sub_D] = size(prop_matrix);

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

save(['./Output/2B' num2str(j) '.mat'], 'store_B', '-v7.3')
save(['./Output/2Theta' num2str(j) '.mat'], 'store_Theta', '-v7.3')

end