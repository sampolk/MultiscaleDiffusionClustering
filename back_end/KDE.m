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
    parfor l=1:n
        p(l)=sum(exp(-(D(1:NN+1,l).^2)./(sigma0^2)));
    end
    
else
    parfor l=1:n
        p(l)=sum(exp(-(D(1:n,l).^2)./(sigma0^2)));
    end
end
p=p./sum(p);
    
end