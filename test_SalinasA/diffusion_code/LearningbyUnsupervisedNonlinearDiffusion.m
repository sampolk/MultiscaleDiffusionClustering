function [C, K, Dt] = LearningbyUnsupervisedNonlinearDiffusion(X, t, G, p, K_known)
%{
 - This function produces a structure with multiscale clusterings produced
   with the LUND algorithm, presented in the following paper. 

        - Maggioni, M., J.M. Murphy. Learning by Unsupervised Nonlinear 
          Diffusion. Journal of Machine Learning Research, 20(160), 
          pp. 1-56. 2019.
    
   and analyzed further in the following papers:

        - Murphy, James M and Polk, Sam L., 2020. A Multiscale Environment 
          for Learning By Diffusion.(In Preparation).
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

Output: 
            - C:                n x 1 vector storing the LUND clustering 
                                of X at time t.
            - K:                Scalar, number of clusters in C.
                                clusters in the Labels(:,t) clustering.
            - Dt:               n x 1 matrix storing \mathcal{D}_t(x). 

© 2021 Sam L Polk, Tufts University. 
email: samuel.polk@tufts.edu
%} 

n = length(X);
C = zeros(n,1);

% Calculate diffusion map
DiffusionMap = zeros(size(G.EigenVecs));
for l = 1:size(DiffusionMap,2)
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
if nargin == 5
    K = K_known;
else
    [~, K] = max(Dt(m_sorting(1:n-1))./Dt(m_sorting(2:n)));
end

if K == 1
    C = ones(n,1);
else

    % Label modes
    C(m_sorting(1:K)) = 1:K;

    % Label non-modal points according to the label of their Dt-nearest
    % neighbor of higher density that is already labeled.
    [~,l_sorting] = sort(p,'descend');
    for j = 1:n
        i = l_sorting(j);
        if C(i)==0 % unlabeled point
            candidates = find(and(p>=p(i), C>0)); % Labeled points of higher density. cannot include xi. 
            [~,temp_idx] = min(DiffusionDistance(i, candidates));
            C(i) = C(candidates(temp_idx));
        end
    end
end
end 