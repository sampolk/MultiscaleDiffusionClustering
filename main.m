%{
Remarks:

 -  You must have run india_forest.m first and downloaded the file titled
    'normalized_data_NNs.mat'. Otherwise, this script will not work. 

 -  On my (Sam's) computer, I can only store weight matrices with 1000 
    nearest neighbors in my workspace due ot memory constraints. My 
    computer has 8 GB of RAM.

%}
%% Load data and nearest neighbors

load('normalized_data_NNs.mat')
Dist_NN = [D1, D2, D3, D4];
Idx_NN = [Idx1, Idx2, Idx3, Idx4];
clear 'D1' 'D2' 'D3' 'D4' 'Idx1' 'Idx2' 'Idx3' 'Idx4'

%% Set Hyperparameters
% 
Hyperparameters.DensityNN = 200; % must be ≤ 3500
Hyperparameters.Sigma0 = prctile(Dist_NN(:,1:Hyperparameters.DensityNN),50, 'all');
Hyperparameters.SpatialParams.ImageSize = [500,500];
Hyperparameters.DiffusionNN = 100;
Hyperparameters.NEigs = 10;
Hyperparameters.DiffusionTime = 0;
Hyperparameters.NumDtNeighbors = 400;

% Endmember algorithm specifications
Hyperparameters.EndmemberParams.Algorithm = 'N-FINDR';
Hyperparameters.EndmemberParams.K = 7; 

Hyperparameters.K_Known = 5; % Optional parameter. If included, it is assigned as the number of clusters. Else, K is found fully unsupervised.

save_on = 0;
plot_on = 1;

%% Calculate Density

p = KDE_large(Dist_NN, Hyperparameters);

if save_on
    save('density.mat','p', 'Hyperparameters')
end
if plot_on
  
    imagesc(reshape(log10(p), M,N))
    colorbar
    xticks([])
    yticks([])
    title('$\log_{10}(p)$', 'interpreter' , 'latex')
    caxis([prctile(log10(p),1), prctile(log10(p),99)])
end

%% Calculate Graph Structure

G = extract_graph_large(X, Hyperparameters, Idx_NN, Dist_NN);

if save_on
    save(strcat('graph_adjesency-', num2str(Hyperparameters.DiffusionNN), '.mat'), 'G')
end

if plot_on 
    for i = 1:min([Hyperparameters.NEigs-1,12])
        subplot(3,4,i)
        imagesc(reshape(G.EigenVecs(:,1+i),500,500))
        pbaspect([1,1,1])
        title(strcat('Eigenvector-', num2str(i+1)))
        xticks([])
        yticks([])
        colorbar
    end
end

%% Calculate LUND Clustering

[C, K, Dt] = LearningbyUnsupervisedNonlinearDiffusion_large(X, Hyperparameters, G, p);

if save_on
    save(strcat('LUND.mat'), 'C', 'K','Dt', 'Hyperparameters')
end
if plot_on 
    imagesc(reshape(C,500,500))
    pbaspect([1,1,1])
    title(strcat('LUND Clustering, $K$=', num2str(K), ', $t$=', num2str(Hyperparameters.DiffusionTime)), 'interpreter', 'latex')
    xticks([])
    yticks([])
    colorbar
end

%% Calculate LUND Clustering with Endmembers

Clustering = LUNDEndmember(X, Hyperparameters, Idx_NN, Dist_NN, G);

if save_on
    save(strcat('LUND_Endmember.mat'), 'Clustering')
end

if plot_on 
    imagesc(reshape(Clustering.Labels,500,500))
    pbaspect([1,1,1])
    title(strcat('LUND Clustering, $K$=', num2str(Clustering.K), ', $t$=', num2str(Hyperparameters.DiffusionTime)), 'interpreter', 'latex')
    xticks([])
    yticks([])
    colorbar
end
