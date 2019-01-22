% Usage: M-step, maximizes the parameters in LDA model which are Dirichlet 
% parameter on both topic distribution for each document (alpha) and 
% word distribution for each topic. (See Blei's paper)

function [alpha,beta]=lda_Mstep(V,alpha_last,phi,gamma_dirc,Bag,lap)

[K,M]=size(gamma_dirc);

beta=zeros(K,V);

words=cell2mat(Bag); % index and counts
clear Bag

% beta
for i=1:size(phi,2)
   curr=words(1,i);
   beta(:,curr)=beta(:,curr)+phi(:,i)*words(2,i) ;   
end

beta=beta+lap; % Smooting ??
beta=normalize(beta,2);


% alpha
alpha_t=alpha_last;
epsilon=0.001;
time=500;

t=0;
e=100;
psiGama=psi(gamma_dirc);
psiSumGama=psi(sum(gamma_dirc,1));

while e>epsilon&&t<time    
    g=sum((psiGama-ones(K,1)*psiSumGama),2)+M*(psi(sum(alpha_t))-psi(alpha_t));    
    h=-M*psi(1,alpha_t);
    z=M*psi(1,sum(alpha_t));
    c=sum(g./h)/(1/z+sum(1./h));
    delta=(g-c)./h;
    
    tao=1;
    alpha_tt=alpha_t-delta;
    while any(alpha_tt<=0)
        tao=tao/2;
        alpha_tt=alpha_t-tao*delta;
    end
    e=sum(abs(alpha_tt-alpha_t))/sum(alpha_t);
    
    alpha_t=alpha_tt;
    
    t=t+1;
end
alpha=alpha_t;