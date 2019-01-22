% Usage: E-step estimates the variational parameters that minimizes the
% KL-diveregence. (See Blei's paper)

function [gamma_dirc phi]=lda_Estep(Bag,Ibag,alpha,beta,phi_old)

M=numel(Ibag);
K=size(alpha,1);
MN=Ibag{end}(end);

gamma_dirc=zeros(K,M);
phi=zeros(K,MN);

if nargin<6
    phi_old=rand(K,MN);
    phi_old=normalize(phi_old,1);
end


for i=1:M
    [gamma_dirc(:,i) phi(:,Ibag{i})]=E_step(Bag{i},alpha,beta,phi_old(:,Ibag{i}));
end

end


function [gamma_dirc phi]=E_step(Bag,alpha,beta,phi_old)

K=numel(alpha);
N=size(phi_old,2);


phi=zeros(K,N);

gamma_dirc_old=ones(K,1);

Iter=500;
i=1;
epsilon=1e-4;
e=100;

ind_nzero=Bag(1,:);
count=Bag(2,:);


%%

while i<Iter && e>epsilon
    
    % gamma_dirc
    gamma_dirc=alpha+phi_old*count'; % Sum along colums
    
    %phi    
    psi_mat=repmat(psi(gamma_dirc)-psi(ones(1,K)*gamma_dirc),1,N);
    
    phi=exp(psi_mat).*beta(:,ind_nzero);
    phi=normalize(phi,1);
    
    %
    e1=diff_vec(gamma_dirc_old,gamma_dirc);
    %e2=diff_vec(phi_old,phi);
    e=e1;
        
    phi_old=phi;
    gamma_dirc_old=gamma_dirc;
    
    i=i+1;
    
end


end