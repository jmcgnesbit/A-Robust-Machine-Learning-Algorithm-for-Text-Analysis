
% Author : F. Yahaya
% Date: 06/09/2018
% Contact: farouk.yahaya@univ-littoral.fr

% clear all variables
% First, initialize the random number generator to make the results in each
% test repeatable.
% using a seed of 1 rng(1)
% Goal: Generate synthetic data of sizes: [m,n]=500, 5000,10000.
% For each of the aforementioned matrix size, a total of 40 tests are made
% with noise of approximately 30db.


% <Parameters>
%       Total_Tests: Number of desired tests.
%       X       : Input data matrix (m x n)
%       r       : Target low-rank
%       nu      : Security rank
%      mat_size : Give a name to your matrix size for easier identification of output

%       Wtheo, Htheo        : Simulation of theoretical matrix Wtheo,Htheo
%       Vtheo               : simulation of theoretical data matrix Vtheo
%       N                   : simulating noise matrice N
%       X = Xtheo+N         : simulating data matrix X
%       Winit,Hinit         : Matrix initialisation
%       SNR                 : Signal to Noise Ratio
%       compressionLevel    : Compression Level, Default=20. See: Mariano
%                             Tepper and Guillermo Sapiro, Compressed Nonnegative
%                             Matrix Factorization is Fast and Accurate, 2015.

clear
clc
rng(1)
                       
Total_Tests =40;

for i =1:Total_Tests
    mat_size='500x500'; % this is just for naming. To change matrix size, change the values for m and n accordingly.                               

    m = 5000; n = 5000; 
    r = 15;
    nu=10;
    
    rng(i)
  
    Wtheo = 10*rand(m,r);
    Htheo = 10*rand(r,n);
    
    Vtheo = Wtheo*Htheo;
    
    rng(200+i)
    Winit = rand(m,r);
    Hinit = rand(r,n);
  
    compressionLevel = 20;
    
    rng(300+i)
    N = 12*randn(m,n); % To reproduce similar results, please keep SNR close to 30db
    SNR = snr(Vtheo,N);
    X=Vtheo+N;
   
    save( ['synthetic_data/data_',mat_size,'_',int2str(i),'.mat'] ,'X','SNR','N' ,'Wtheo' , 'Htheo' , 'Winit' , 'Hinit' , 'Vtheo','compressionLevel','nu','r' )
    
    
end
disp(['matrix size ', '( ',mat_size,' )', ' : all ',int2str(Total_Tests), ' data simulations ', ' completed!' ]);


