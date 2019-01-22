% Usage: The main function to learn LDA model on your data. Based on
% variational formulation in Blei's Original paper. 

% Input: 
% 1. X: Sparse matrix representation BoW counts in N*D format, where N is
% no of points and D is no of words.
%2. K: no of topics
%3. Words: cell containing all the words. Use to show the topics learned.

%Output: 
%1. alpha: Dirichlet-parameter for per-document topic distribution (No_topics x 1 matrix)
%2. Beta: Learned per-topic word distribution (No_topics x words matrix) 
%3. Phi: Variational parameter (for Z_{mn}) denoting topic-wise distribution for
% all words appearing in all documents (No_topics x enitre_words matrix) 
%4. gamma_dirc: topic-wise distribution for each document (variational
%parameter) (No_topics*No_of_documents matrix)


%Comment: Also outputs the top words for each learned topic in different iterations
% in the folder ./topics
% Pl. see README for a description.

% Karan Sikka, UCSD, 8/10/2013
% contact: ksikka@ucsd.edu/karan.sikka1@gmail.com,http://mplab.ucsd.edu/~ksikka/

function [alpha,beta,phi,gamma_dirc]=learn_lda(X,K,words)


%% Initial Parameters


epsilon=0.00001; % Stopping criteria 
max_iter=500; % Adjust iterations
lap=0.00001; % Hyperprior for Beta

%% Input and Output Setting
% Create topics folder again
if ~exist('topics','dir') % Output top-words in each topic
    mkdir('topics');
else
    cd ('topics')
    delete '*.txt'
    cd ..
end

M=size(X,1);
fprintf('Topics = %d\n',K)
V=size(X,2);


% Converting X in a particular format (see below)
Bag=cell(1,M); % Stores indices and counts of non-zero words in each document seperately
Ibag=cell(1,M); % Indicize each word in a document to build phi matrix
k=1;

for i=1:M
    ind_nzero=find(X(i,:));
    Bag{i}=full([ind_nzero;X(i,ind_nzero)]);
    Ibag{i}=k:k+numel(ind_nzero)-1;
    k=k+numel(ind_nzero);
end

clear X

MN=k-1; % non-zero words
clear k

%% Initializing variables
alpha_old=rand(K,1);
beta_old=rand(K,V);
beta_old=normalize(beta_old,2);
phi_old=single(ones(K,MN)/K);


disp('learning LDA')
lik=0; % Stores likelihood

%%
for glb=1:max_iter
    
    
    % E-step
    fprintf('%d. \nE Step',glb)
    tic
    [gamma_dirc,phi]=lda_Estep(Bag,Ibag,alpha_old,beta_old,phi_old); % Expectation-Step
    t=toc;
    fprintf('............... %.1f sec\n',t);
    
    
    % M-step
    fprintf('M Step')
    tic
    [alpha,beta]=lda_Mstep(V,alpha_old,phi,gamma_dirc,Bag,lap); % Maximization step
    t=toc;
    fprintf('............... %.1f sec\n',t);
    
    
    alpha_old=alpha; % Updating
    beta_old=beta;
    
    
    lik(glb)=lda_eval_lik(Bag,Ibag,phi,beta,gamma_dirc,alpha); % Evaluating Likelihood
    fprintf('Lik   ............... %.2f\n\n\n',lik(glb));
    
    if glb>1 && converged(lik(glb),lik(glb-1),epsilon) % Check for convergence
        break
    end
    
    for i=1:K
        f=fopen(['./topics/topics_' num2str(i) '.txt'],'a');
        [a b]=sort(beta(i,:),'descend');
        fprintf(f,'%d. %s ',glb,words{b(1)});
        for j=2:30
            fprintf(f,'%s ',words{b(j)});
        end
        fprintf(f,'\n \n \n');
        fclose(f);
    end
    
    
end

end


