function G = extract_graph_large(X, Hyperparameters, Idx_NN, Dist_NN)
%{
 - This function produces a graph structure to be used in cluster analysis
   on large graphs. 

 - If spatial parameters are provided, this function automatically computes
   a modified graph that directly incorporates spatial geometry into graph
   diffusion. See paper below for more. If spatial parameters aren't
   provided, this function reverts to a standard KNN graph.

        - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-
          Spatial Diffusion Geometry for Hyperspectral Image Clustering. 
          (In Review)

- If the NEigs is not included in the Hyperparameters structure, m = argmax( abs(diff(EigenVals)))
  eigenvalues are included in the Graph structure.

Inputs: X:                      (M*N)xD Data matrix .
        Hyperparameters:        Structure with graph parameters with the 
                                required fields: 
            - SpatialParams:    Stores the dimensions of the original image
            - DiffusionNN:      Number of nearest neighbors in KNN graph.
            - WeightType:       Equal to either 'adjesency' or 'gaussian' (Optional).
            - Sigma:            If WeightType == 'gaussian', then diffusion scale parameter Sigma>0 required.

        Idx_NN:             Indices of l2-nearest neighbors of points.
        Dist_NN:            Distances between points and their l2-nearest neighbors.

Output: Graph structure with the following fields:
            - Hyperparameters:  Parameters used to generate Graph structure.
            - EigenVecs:        First m eigenvectors of P.
            - EigenVals:        First m eigenvalues of P.

© 2021 Sam L Polk, Tufts University. 
email: samuel.polk@tufts.edu
%}

n = length(X); 
NN = Hyperparameters.DiffusionNN;
n_eigs = Hyperparameters.NEigs;


if ~isfield(Hyperparameters, 'WeightType')
    Hyperparameters.WeightType = 'adjesency';
end
if strcmp(Hyperparameters.WeightType, 'gaussian')
    Hyperparameters.Sigma = median(Dist_NN(:,1:NN),'all');
end

if isfield(Hyperparameters.SpatialParams, 'SpatialRadius') 

    W = spatial_weight_matrix(X,Hyperparameters);
    
else
    
    % Preallocate memory for sparse matrix computation
    ind_row = repmat((1:n)', 1,NN);  % row indices for nearest neighbors.
    ind_col = Idx_NN(:,1:NN);         % column indices for nearest neighbors.
    if strcmp(Hyperparameters.WeightType, 'adjesency')

        W = sparse(ind_row, ind_col, ones(n,NN)); % construct W
        W = adjacency((graph((W+W')./2)));   % Convert W to adjesency matrix. (W+W)./2 forces symmetry.
    elseif strcmp(Hyperparameters.WeightType, 'gaussian')

        if ~isfield(Hyperparameters, 'Sigma')
            sigma = prctile(Dist_NN(:,1:NN),50,'all');
        else
            sigma = Hyperparmeters.Sigma;
        end

        W = sparse(ind_row, ind_col, exp(-(Dist_NN(:,1:NN).^2)./(sigma^2))); % construct W
        W = (W+W')./2;   %  (W+W)./2 forces symmetry
    end
end

Dinv = spdiags(1./sum(W)',0,n,n);    % D^{-1}.

[V,D, flag] = eigs(Dinv*W, n_eigs, 'largestabs'); % First n_eigs eigenpairs of transition matrix

 
if flag
    disp('Convergence Failed.')
    G = NaN;
else
    [lambda,idx] = sort(diag(abs(D)),'descend');
    lambda(1) = 1;
    V = real(V(:,idx));

    G.Hyperparameters = Hyperparameters;
    G.EigenVecs = V;
    G.EigenVals = lambda;
end

