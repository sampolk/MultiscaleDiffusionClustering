function [p,D] = KDE(X,Hyperparameters, D)

% Computes the local densitys of the points in an n x D matrix X, using the
% threshold value th.  Use K nearest neighbor.  A matrix of distances may
% be passed in.  

% Extract hyperparameters
NN = Hyperparameters.DensityNN;
sigma0 = Hyperparameters.Sigma0;

n = length(X);

% Initialize
p=zeros(n,1);

if nargin<3
    D=pdist(X);
    D=squareform(D); %Pairwise distances between points in X
end

D=sort(D);

if NN<n
    p = sum(exp(-(D(1:NN+1,:).^2)./(sigma0^2)))';
else
    p = sum(exp(-(D.^2)./(sigma0^2)))';
end
p=p./sum(p);
    
end