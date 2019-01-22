function lik=lda_eval_lik(Bag,Ibag,phi,beta,gamma_dirc,alpha) % Evaluates likelihood

M=numel(Ibag);
K=size(phi,1);

words=cell2mat(Bag);
clear Bag;

% Precompute gamma_hat
gamma_hat=ones(1,K)*gamma_dirc;

% Term 1
t1=-ones(1,K)*(phi.*log(phi+eps))*words(2,:)';


% Term 2
t2=0;
for i=1:M
    indx_elem=Ibag{i};    
    fact_psi=psi(gamma_dirc(:,i))-psi(gamma_hat(i)); % The psi factor for each bag
    t2=t2+ones(1,K)*(phi(:,indx_elem).*repmat(fact_psi,[1 numel(Ibag{i})]))*words(2,indx_elem)';
end
clear indx_elem fact_psi

% Term 3
t3=ones(1,K)*(phi.*log(beta(:,words(1,:))))*words(2,:)';
clear phi


% Term4
t4=0;
for i=1:M    
    t4=t4+ones(1,K)*gammaln(gamma_dirc(:,i))-gammaln(gamma_hat(i))-(gamma_dirc(:,i)-1)'*(psi(gamma_dirc(:,i))-psi(gamma_hat(i)));
end

% Term 5
t5=0;
for i=1:M        
    d2=gammaln(ones(1,K)*alpha)-ones(1,K)*gammaln(alpha);        
    d1= (alpha-1)'*(psi(gamma_dirc(:,i))-psi(gamma_hat(i)));
    t5=t5+d1+d2;    
end



lik=(t1+t2+t3+t4+t5);
t1;
t2;
t3;
t4;
t5;

end