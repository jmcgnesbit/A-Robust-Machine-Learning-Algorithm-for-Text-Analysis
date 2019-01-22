C= doc2cell(cleanDocuments);

wordvec = string([]);

for d = 1:D
    wordvec = [wordvec; C{d}'];
end
wordint = grpcategorical(wordvec)
bow.Vocabulary;
tdmat
N_K = 40;


alpha = 1;
gamma = 1;

% Word-topic assignment
z = zeros(D, maxwords); 
for i = 1:N_D
    for l = 1:N_W(i)
        z(i,l) = randi(N_K);
    end
end

% Document-topic assignment
Pi = zeros(N_D,N_K);

for i = 1 : D
    Pi(i, :) = drchrnd(alpha*ones(K,1)', 1);
end

B = zeros(N_K, N_W);

for k = 1:N_K
    B(k, :) = drchrnd(gamma*ones(K,1)', 1);
end
