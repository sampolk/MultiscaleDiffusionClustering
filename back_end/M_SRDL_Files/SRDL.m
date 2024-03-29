function [C, K, Dt] = SRDL(X, t, Hyperparameters, G, p)
%{
 - This function produces a structure with multiscale clusterings produced
   with the SRDL algorithm, presented in the following paper. 

        - Murphy, James M., and Mauro Maggioni. "Spectral-spatial diffusion 
          geometry for hyperspectral image clustering." IEEE Geoscience and 
          Remote Sensing Letters (2019).

   and analyzed further in the following paper:

        - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-
          Spatial Diffusion Geometry for Hyperspectral Image Clustering. 
          (In Review)

Inputs: X:                      Data matrix.
        t:                      Diffusion Time Step
        Hyperparameters:        Optional structure with graph parameters
                                with required fields:  
                                    - Sigma>0: diffusion scale.
                                    - DiffusionNN: no. nearest neighbors in
                                                   KNN graph.
                                    - NEigs: no.   eigenvectors. (Optional)
                                    - Sigma0>0:    KDE Bandwidth.
                                    - DesnityNN:   no. nearest neighbors  
                                                   used in KDE.
        G:                      Graph structure computed using  
                                'extract_graph.m' (Optional).
        p:                      Kernel Density Estimator (Optional).

Output: Clusterings structure with the following fields:

            - C:                n x 1 vector storing the SRDL clustering 
                                of X at time t.
            - K:                Scalar, number of clusters in C.
                                clusters in the Labels(:,t) clustering.
            - Dt:               n x 1 matrix storing \mathcal{D}_t(x). 

© 2021 Sam L Polk, Tufts University. 
email: samuel.polk@tufts.edu
%}


if nargin == 3
    G = extract_graph(X, Hyperparameters);
    p = KDE(X, Hyperparameters);
elseif nargin == 4
    p = KDE(X, Hyperparameters);
end    

n = length(X);
C = zeros(n,1);

% Calculate diffusion map

if isfield(Hyperparameters, 'NEigs')
    n_eigs = Hyperparameters.NEigs;
else
   [~, n_eigs] = max(abs(diff(G.EigenVals)));
end

DiffusionMap = zeros(n,n_eigs);
for l = 1:n_eigs
    DiffusionMap(:,l) = G.EigenVecs(:,l).*(G.EigenVals(l).^t);
end

% Calculate pairwise diffusion distance at time t between points in X
DiffusionDistance = squareform(pdist(real(DiffusionMap)));

% compute rho_t(x), stored as rt
rt = zeros(n,1);
for i=1:n
    if ~(p(i)==max(p))
        rt(i)=min(DiffusionDistance(p>p(i),i));     
    else
        rt(i)=max(DiffusionDistance(i,:)); 
    end
end

% Extract Dt(x) and sort in descending order
Dt = rt.*p;
[~, m_sorting] = sort(Dt,'descend');

% Determine K based on the ratio of sorted Dt(x_{m_k}). 
if isfield(Hyperparameters, 'K_known')
    K = K_known;
else
    [~, K] = max(Dt(m_sorting(1:n-10))./Dt(m_sorting(2:n-9))); % avoid some trivial clusterings by deleting the last 9 entries of the sorted values
end

if K == 1
    C = ones(n,1);
else
    
    % Label modes
    C(m_sorting(1:K)) = 1:K;

    % Label non-modal points according to the label of their Dt-nearest
    % neighbor of higher density that is already labeled.
    [~,idx] = sort(p,'descend');
    
    % Define spatial consensus variables
    M = Hyperparameters.SpatialParams.ImageSize(1);
    N = Hyperparameters.SpatialParams.ImageSize(2);  
    WindowSize = Hyperparameters.SpatialParams.ConsensusSpatialRadius;  

    if WindowSize == 0   
        % Equivalent to labeling pixels according to their Dt-nearest
        % neighbor of higher density that is already labeled
        for j = 1:n
            i = idx(j);
            if C(i)==0 % unlabeled point
                candidates = find(and(p>=p(i), C>0)); % Labeled points of higher density.
                [~,temp_idx] = min(DiffusionDistance(i, candidates));
                C(i) = C(candidates(temp_idx));
            end
        end
    else
        %We incorporate spatial consensus in this case. 
        [I,J]=find(sum(reshape(X,M,N,[]),3)~=0); % indices of each pixel
    
        % First pass
        for j=1:length(C)
            if C(idx(j))==0
    
                % Find diffusion distance nearest neighbor of higher density
                % that is already labeled
                [~,NN]=min(DiffusionDistance(idx(j),idx(1:j-1)));
                C(idx(j))=C(idx(NN));
    
                % Construct spatial consensus label
                neighbors = FindNeighbors([I(idx(j)),J(idx(j))],WindowSize,M,N)'; % indices of spatial neighbors of X(idx(j),:).
                neighboringPixels = sub2ind([M,N], neighbors(:,1), neighbors(:,2)); % Convert to indices 
                labeledNeighboringPixels = intersect(neighboringPixels, idx(1:j-1)); % contains the indices of neighboring pixels that are already labeled
    
                if ~isempty(labeledNeighboringPixels) % True if no spatial consensus label. 
                    SpatialConsensusLabel = mode(C(labeledNeighboringPixels));
                    if ~(SpatialConsensusLabel  == C(idx(j)))
                        C(idx(j)) = 0;
                    end
                end
            end
        end
            
        % Second pass
        for j=1:length(C)
            if C(idx(j))==0
                % Find diffusion distance nearest neighbor of higher density
                % that is already labeled
                [~,NN]=min(DiffusionDistance(idx(j),idx(1:j-1)));
                C(idx(j))=C(idx(NN)); 
    
                % Construct spatial consensus label
                neighbors = FindNeighbors([I(idx(j)),J(idx(j))],WindowSize,M,N)'; % indices of spatial neighbors of X(idx(j),:).
                neighboringPixels = sub2ind([M,N], neighbors(:,1), neighbors(:,2)); % Convert to indices 
                labeledNeighboringPixels = intersect(neighboringPixels, find(C>0)); % contains the indices of neighboring pixels that are already labeled
    
                if ~isempty(labeledNeighboringPixels) % True if a spatial consensus label exists for X(idx(j),:).
                    [SpatialConsensusLabel, SpatialConsensusFrequency] = mode(C(labeledNeighboringPixels));  
                    if SpatialConsensusFrequency>0.5*length(labeledNeighboringPixels) % true if the majority of labeled pixels neighboring X(idx(j),:) have label = SpatialConsensusFrequency
                        C(idx(j))=SpatialConsensusLabel;
                    end 
                end
            end
        end
    end
end