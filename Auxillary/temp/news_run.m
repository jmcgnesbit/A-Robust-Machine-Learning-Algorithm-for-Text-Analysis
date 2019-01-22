load('./20news-bydate/20news_train_bow.mat')
load('./20news-bydate/20news_words.mat')
[alpha,beta,phi,gamma_dirc]=learn_lda(X,110,words)