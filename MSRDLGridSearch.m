%% 

NNs = [10:10:100];
prctiles = 5:10:95;
Rs = 2:15;

%% 
[X,Y] = extract_salinasA(); %load data
load('salinasA-SRDL-HP.mat')
Hyperparameters = rmfield(Hyperparameters,'Sigma');

X = X./vecnorm(X,2,2);
X = X + (10^(-7)).*randn(size(X));

% Calculate 1000 nearest neighbors
[Idx_NN, Dist_NN] = knnsearch(X,X,'K', 1001);
Idx_NN(:,1) = [];
Dist_NN(:,1) = [];

Dists = squareform(pdist(X));
DistsSorted = sort(Dists);

%% Run grid search

Ks = cell(length(NNs), length(Rs), length(Rs), length(prctiles));

k=0;
l=0;
for i = 1:length(NNs)
    for j = 1:length(Rs)

        z = zeros(M*N,1);
        for l = 1:M*N
            [idx1,idx2] = ind2sub([M,N], l);
            z(l) = size(FindNeighbors([idx1,idx2], Rs(j), M, N),2)-1; % number of spatial neighbors at each pixel, minus 1 to exclude X(l,:).
        end
        l=0;

        if NNs(i)<=min(z)

            for k = 1:length(Rs)

                for l = 1:length(prctiles)

                    % Set relevant hyperparameters
                    Hyperparameters.DiffusionNN = NNs(i);
                    Hyperparameters.DensityNN = NNs(i);
                    Hyperparameters.SpatialParams.GraphSpatialRadius=Rs(j);
                    Hyperparameters.SpatialParams.ConsensusSpatialRadius=Rs(k);

                    % Extract graph
                    G_SpatialRegularization = extract_graph(X, Hyperparameters);
            
                    if ~isnan(G_SpatialRegularization.EigenVals(1)) && G_SpatialRegularization.EigenVals(2)<1 && min(G_SpatialRegularization.StationaryDist) >0
                        % True if eigendecomposition converges to reasonable graph
        
                        Hyperparameters.Sigma0 = prctile(Dist_NN,prctiles(l),'all');
                        p = KDE(X,Hyperparameters, DistsSorted);
        
                        if sum(isnan(p)) == 0 && min(p)>0
            
                            Cs = M_SRDL(X, Hyperparameters, G_SpatialRegularization, p);
                            Ks{i,j,k,l} = Cs.K;

                        end
                    end
                    disp([i,j,k,l]./size(nmis))        
                end
            end
        end
        disp([i,j,k,l]./size(nmis))        
    end
end



        

        