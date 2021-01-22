function C = SpectralClustering(K,G)

n = length(G.EigenVecs);

U = G.EigenVecs(:,1:K);
Y = zeros(n,K);

for i = 1:n
    Y(i,:) = U(i,:)./norm(U(i,:),2);
end
C = kmeans(Y, K);