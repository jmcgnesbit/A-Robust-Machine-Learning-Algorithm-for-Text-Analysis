% Author : F. Yahaya
% Date: 06/09/2018
% Contact: farouk.yahaya@univ-littoral.fr

clear
clc
Total_Tests=40;
mat_size='500x500'; % this is just for naming. To change matrix size, open the "data_simulation.m" file and change                                
                        % your mxn matrix size as desired.

for i=1:Total_Tests
    
  load(['output/RSI_NeNMF_',mat_size,'_',int2str(i),'.mat'] , 'RRE_RSI_NeNMF', 'T_RSI_NeNMF')
%     load(['output/RPI_',mat_size,'_',int2str(i),'.mat'] , 'RRE_RPI','T_RPI' )
    load(['output/vanilla_NeNMF_',mat_size,'_',int2str(i),'.mat'] , 'RRE_VANILLA_NeNMF', 'T_VANILLA_NeNMF' )
    
    
    rRE_RSI_NeNMF(i,:)  = RRE_RSI_NeNMF;
    t_RSI_NeNMF(i,:) = T_RSI_NeNMF;
    
%     rRE_RPI(i,:)  = RRE_RPI;
%     t_RPI(i,:) = T_RPI;
    rRE_VANILLA_NeNMF(i,:)  = RRE_VANILLA_NeNMF;
    t_VANILLA_NeNMF(i,:) = T_VANILLA_NeNMF;
    
end

figure,
subplot(211)
semilogy(t_VANILLA_NeNMF',rRE_VANILLA_NeNMF','b')
hold on, semilogy(t_RSI_NeNMF',rRE_RSI_NeNMF','r')
axis([0 15 1e-4 1.01]) 
xlabel('CPU time (s)')
ylabel('RRE')




subplot(212), semilogy([0:length(t_VANILLA_NeNMF)-1],rRE_VANILLA_NeNMF','b')
hold on, semilogy([0:length(t_RSI_NeNMF)-1],rRE_RSI_NeNMF','r')
axis([0 40 1e-4 1.01]) 
xlabel('Iterations')
ylabel('RRE')

