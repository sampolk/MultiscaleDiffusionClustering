function [X,Y] = extract_salinasA()

load('SalinasA_smallNoise.mat');
load('SalinasA_gt.mat');

[M,N] = size(salinasA_gt);
Y=double(salinasA_gt);
X=permute(X,[2,1,3]);


X=reshape(X,size(X,1)*size(X,2),size(X,3));
X=X./repmat(sqrt(sum(X.*X,1)),size(X,1),1); %Normalize X

n = length(X); 

Y_rshp = reshape(Y,M*N,1);

GT = zeros(n,1);
unique_k = unique(Y_rshp);
K_GT = length(unique_k);
for k = 1:K_GT
    GT(Y_rshp==unique_k(k)) = k;
end
Y = GT;