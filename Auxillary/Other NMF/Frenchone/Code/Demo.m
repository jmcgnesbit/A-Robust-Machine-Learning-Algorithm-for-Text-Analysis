% Author : F. Yahaya
% Date: 06/09/2018
% Contact: farouk.yahaya@univ-littoral.fr

% clear all variables
% First, initialize the random number generator to make the results in each
% test repeatable.
% using a seed of 1 rng(1)

% Goal: Execute the algorithms viz., Vanila NeNMF, Randomized Power
% Iterations NeNMF (RPI), Randomized Subspace Iterations NeNMF (RSI) for
% sizes: [m,n]=500, 5000,10000.

% Tmax : Execution time in seconds.
% mat_size : Give a name to your matrix size for easier identification of output
% Total_Tests : Total number of test to perform. E.g. we made 40 tests.

clear
clc
Tmax =15;
Total_Tests=40;



for i=1:Total_Tests
    
   mat_size='500x500'; % this is just for naming. To change matrix size, open the "data_simulation.m" file and change                                
                        % your mxn matrix size as desired.
    
    
%     % Randomized Subspace Iterations (NeNMF) 
    
    clearvars -except i Tmax mat_size Total_Tests;
   load( ['synthetic_data/data_',mat_size,'_',int2str(i),'.mat'])
    [ W_RSI_NeNMF , H_RSI_NeNMF,RRE_RSI_NeNMF, T_RSI_NeNMF] = RSI_NeNMF(X, Winit , Hinit ,r, Tmax);
   save( ['output/RSI_NeNMF_',mat_size,'_',int2str(i),'.mat'], 'RRE_RSI_NeNMF', 'T_RSI_NeNMF', '-v7.3' )
  
    % VANILA NeNMF
    clearvars -except i Tmax mat_size Total_Tests;
   load( ['synthetic_data/data_',mat_size,'_',int2str(i),'.mat'])
    [ W_VANILLA , H_VANILLA,RRE_VANILLA_NeNMF, T_VANILLA_NeNMF ] = VANILLA_NeNMF(X , Winit , Hinit,Tmax );
   save( ['output/vanilla_NeNMF_',mat_size,'_',int2str(i),'.mat'], 'RRE_VANILLA_NeNMF', 'T_VANILLA_NeNMF', '-v7.3' )

    disp( ['iter:',int2str(i), ' of ', mat_size , '   OK!'] )
    
end

disp(['matrix size ', '( ',mat_size,' )', ' : all ',int2str(Total_Tests), ' Tests', ' completed!' ]);


