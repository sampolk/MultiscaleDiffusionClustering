function [C, K, Dt] = LearningbyUnsupervisedNonlinearDiffusion_large(X, Hyperparameters, G, p)
%{
 - This function produces a structure with multiscale clusterings produced
   with the LUND algorithm, presented in the following paper. 

        - Maggioni, M., J.M. Murphy. Learning by Unsupervised Nonlinear 
          Diffusion. Journal of Machine Learning Research, 20(160), 
          pp. 1-56. 2019.
    
   and analyzed further in the following papers:

        - Murphy, James M and Polk, Sam L., 2020. A Multiscale Environment 
          for Learning By Diffusion. arXiv preprint, arXiv:2102.00500.
        - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-
          Spatial Diffusion Geometry for Hyperspectral Image Clustering. 
          (In Review)

Inputs: X:                      Data matrix.
        Hyperparameters:        Optional structure with graph parameters
                                with required fields:  
                                    - DiffusionTime:    Diffusion time parameter.
        G:                      Graph structure computed using  
                                'extract_graph_large.m' 
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

if ~isfield(Hyperparameters, 'NEigs')
    Hyperparameters.NEigs = size(G.EigenVecs,2);
end

n = length(X);
C = zeros(n,1);

% Calculate diffusion map
DiffusionMap = zeros(size(G.EigenVecs,1), Hyperparameters.NEigs);
parfor l = 1:size(DiffusionMap,2)
    DiffusionMap(:,l) = G.EigenVecs(:,l).*(G.EigenVals(l).^Hyperparameters.DiffusionTime);
end

disp('Diffusion Map Calculated')

% compute rho_t(x), stored as rt
rt = zeros(n,1);
p_max = max(p);
% p_pct = prctile(p, Hyperparameters.PPct);
parfor i=1:n
    if p(i) == p_max
        rt(i) = max(pdist2(DiffusionMap(i,:),DiffusionMap));
    else
        rt(i) = min(pdist2(DiffusionMap(i,:),DiffusionMap(p>p(i),:)));
    end
        
    if mod(i,100) == 0
        disp(strcat('Rho calculation, ', num2str((1-i/n)*100, 3), '% complete.'))
    end
end

% Extract Dt(x) and sort in descending order
Dt = rt.*p;
[~, m_sorting] = sort(Dt,'descend');


% Determine K based on the ratio of sorted Dt(x_{m_k}). 
if isfield(Hyperparameters, 'K_Known')
    K = Hyperparameters.K_Known;
else
    [~, K] = max(Dt(m_sorting(1:n-1))./Dt(m_sorting(2:n)));
end
disp('Cluster Modes Calculated')

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
            candidates = find(and(p>=p(i), C>0));% Labeled points of higher density. 
            Dtxi = pdist2(DiffusionMap(i,:),DiffusionMap(candidates,:));
            [~,temp_idx] = min(Dtxi);
            C(i) = C(candidates(temp_idx));
        end
        if mod(j,100) == 1
            disp(strcat('Non-modal labeling, ', num2str((j/n)*100, 3), '% complete.'))
        end
        
    end
end 