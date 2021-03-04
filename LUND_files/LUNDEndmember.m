function Clustering = LUNDEndmember(X, Hyperparameters, Idx_NN, Dist_NN, G)
%{
 - This function produces a structure with clustering produced with the 
   LUND algorithm, presented in the following paper. 

        - Maggioni, M., J.M. Murphy. Learning by Unsupervised Nonlinear 
          Diffusion. Journal of Machine Learning Research, 20(160), 
          pp. 1-56. 2019.
    
   and analyzed further in the following papers:

        - Murphy, James M and Polk, Sam L., 2020. A Multiscale Environment 
          for Learning By Diffusion. arXiv preprint, arXiv:2102.00500.
        - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-
          Spatial Diffusion Geometry for Hyperspectral Image Clustering. 
          (In Review)

 - Endmember extraction is incorporated using the algorithm specified in
   Hyperparameters.EndmemberParams.Algorithm. As of now, the number of
   endmembers must be specified as Hyperparameters.EndmemberParams.K. This
   approach is modified from the following paper: 

        - Kangning Cui and Robert J Plemmons. Unsupervised Classication of 
          AVIRIS-NG Hyperspectral Images. In: 2021 WHISPERS (2021).

Inputs:     X:                  (M*N)xD Data matrix .
            Hyperparameters:    Structure with the following fields:
                - SpatialParams:    Stores the dimensions of the original image. 
                - DiffusionTime:    Diffusion time parameter.
                - DiffusionNN:      Number of nearest neighbors in KNN graph.
                - EndmemberParams:  Endmember extraction algorithm and number of endmembers (Optional).
                - WeightType:       Equal to either 'adjesency' or 'gaussian' (Optional).
                - Sigma:            If WeightType == 'gaussian', then diffusion scale parameter Sigma>0 required.
                - IncludeDensity:   1 if density is to be included, 0 otherwise (Optional).
                - DensityNN:        If IncludeDensity == 1, then number of nearest neighbors to comput KDE is required.
                - Sigma0:           If IncludeDensity == 1, then KDE bandwidth Sigma0>0 is required.
            Idx_NN:             Indices of l2-nearest neighbors of points.
            Dist_NN:            Distances between points and their l2-nearest neighbors.
            G:                  Precomputed graph structure (Optional).

Outputs:    Clustering Structure with the following fields:
                - Labels:           LUND+Endmember Clustering of X.
                - K:                Number of clusters.
                - Dt:               \mathcal{D}_t(x).
                - Hyperparameters:  Inputted Hyperparameter structure.
                - Endmembers:       Calculated endmembers.
                - PixelPurity:      Purity of pixels, as measured by endmember unmixing.
%}

if nargin == 4
    % Compute Graph
    G = extract_graph_large(X, Hyperparameters, Idx_NN, Dist_NN);
end

% Compute pixel purity
[p, endmembers] = compute_purity(X,Hyperparameters);

if isfield(Hyperparameters, 'IncludeDensity')
    if Hyperparameters.IncludeDensity

        % Compute Density
        density = KDE_large(Dist_NN, Hyperparameters);
        
        % Store purity and density
        Clustering.PixelPurity = p;
        Clustering.Density = density;

        p = harmmean([p./max(p), density./max(density)],2);
        
    else 
        Clustering.PixelPurity = p;
    end
else
    Clustering.PixelPurity = p;
end

[C, K, Dt] = LearningbyUnsupervisedNonlinearDiffusion_large(X, Hyperparameters, G, p);

% Store Results
Clustering.Labels = C;
Clustering.K = K;
Clustering.Dt = Dt;
Clustering.Hyperparameters = Hyperparameters;
Clustering.Endmembers = endmembers;
Clustering.Graph = G;
