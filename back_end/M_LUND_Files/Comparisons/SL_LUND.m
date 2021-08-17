function [C_Orig, C_MS] = SL_LUND(X, Clusterings, flag)

if nargin == 2
    flag=1;
end

n = length(X);

t_idx = find(and(Clusterings.K>1, Clusterings.K<n/2), 1, 'first');
t = Clusterings.TimeSamples(t_idx);
C_Orig = Clusterings.Labels(:,t_idx); % Single Scale LUND clustering
K = length(unique(C_Orig));

% Compute Pairwise distances
if flag
    
    % Construct Diffusion Distances
    DiffusionMap = zeros(n,Clusterings.Hyperparameters.NEigs);
    for k = 1:size(DiffusionMap,2)
        DiffusionMap(:,k) = ((Clusterings.Graph.EigenVals(k)).^t).*Clusterings.Graph.EigenVecs(:,k);
    end
    Distance = squareform(pdist(DiffusionMap));
    
else
    % Euclidean Distances
    Distance = squareform(pdist(X));    
end

% Implement SLC with using Distances, intiialized at LUND clustering
C_MS  = ones(n, K); 
C_MS(:,1) = C_Orig;
for l = 2:K-1

    % Clustering to merge clusters from
    Cl = C_MS(:,l-1);
    Kl = length(unique(Cl));

    % Compute minimum pairwise diffusion distances
    Distsbtw = NaN*zeros(Kl);
    for k1 = 1:Kl
       for k2 = setdiff(1:Kl, k1)
           Distsbtw(k1,k2) = min(Distance(Cl == k1, Cl == k2), [], 'all');              
       end
    end

    % Find closest two clusters
   [~,i] = min(Distsbtw(:));
   [k1, k2] = ind2sub([Kl, Kl], i);

   % Merge clusters
    Ctemp = Cl;
    Ctemp(Cl == k1) = k2;
    unique_C = unique(Ctemp);
    for k = 1:length(unique_C)
        Ctemp(Ctemp == unique_C(k)) = k;
    end
    C_MS(:,l) = Ctemp;

end
 