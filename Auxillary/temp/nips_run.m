load('nips_bow.mat')
load('nips_words.mat')
[alpha,beta,phi,gamma_dirc]=learn_lda(X,50,words)