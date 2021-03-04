function [C, K, Dt] = SRDL(X, t, Hyperparameters, G, p)

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
    K = Hyperparameters.K_known;
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
    
    % Define spatial consensus variables
    M = Hyperparameters.SpatialParams.ImageSize(1);
    N = Hyperparameters.SpatialParams.ImageSize(2);  
    R = Hyperparameters.SpatialParams.SpatialRadius;  
    [I,J]=find(sum(X,3)~=0); 
    
    % Labeling Pass 1
    for j = 1:n
        i = l_sorting(j);
        if C(i)==0 % unlabeled point
            candidates = find(and(p>=p(i), C>0)); % Labeled points of higher density.
            [~,temp_idx] = min(DiffusionDistance(i, candidates));
            C(i) = C(candidates(temp_idx));
        end
    end
    
    % Labeling Pass 2
    for j=1:n
        
        i = l_sorting(j);
        if C(i)==0
            
            candidates = find(and(p>=p(i), C>0)); % Labeled points of higher density.
            [~,temp_idx] = min(DiffusionDistance(i, candidates)); % index of the Dt-nearest neighbor of higher density that is already labeled. 
            C(i) = C(candidates(temp_idx));
            
            % Special case when multiple points have equal density
            if C(i) == 0
                temp = find(Labels>0); % indices of labeled points
                NN = knnsearch(X(temp,:),X(i,:)); 
                C(i)=C(temp(NN));% Assign label of nearest Euclidean distance-nearest neighbor that is already labeled.
            end
            
            % Check spatial consensus
            try
                LabelsMatrix=sparse(I(l_sorting(1:j)),J(l_sorting(1:j)),C(l_sorting(1:j)),M,N);
                SpatialConsensusLabel=SpatialConsensus(LabelsMatrix,i,M,N,I,J,R);
                if ~(SpatialConsensusLabel==C(i)) && (SpatialConsensusLabel>0)
                    C(l_sorting(j))=0; % reset if spatial consensus label disagrees with DL label
                end
            catch
                keyboard
            end
        end
    end


    % Second pass
    for j=1:n
        i = l_sorting(j);
        
        if C(i)==0
            
            candidates = find(and(p>=p(i), C>0)); % Labeled points of higher density.
            [~,temp_idx] = min(DiffusionDistance(i, candidates)); % index of the Dt-nearest neighbor of higher density that is already labeled. 
            C(i) = C(candidates(temp_idx));
            
            % Special case when multiple points have equal density
            if C(i) == 0
                temp = find(Labels>0); % indices of labeled points
                NN = knnsearch(X(temp,:),X(i,:)); 
                C(i)=C(temp(NN));% Assign label of nearest Euclidean distance-nearest neighbor that is already labeled.
            end
            
            % Check spatial consensus
            try
                LabelsMatrix=sparse(I(l_sorting(1:j)),J(l_sorting(1:j)),C(l_sorting(1:j)),M,N);
                [SpatialConsensusLabel,ConsensusFrac]=SpatialConsensus(LabelsMatrix,i,M,N,I,J,R);
                if ConsensusFrac>.5
                    C(i)=SpatialConsensusLabel;
                end
            catch
                keyboard
            end
        end
    end
    
end
end