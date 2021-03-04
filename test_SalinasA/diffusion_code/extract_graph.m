function Graph = extract_graph(X, Hyperparameters, Dist)
%{
 - This function produces a graph structure to be used in cluster analysis.
   The graph it constructs is a KNN graph with Gaussian kernel edge weights 

 - If spatial parameters are provided, this function automatically computes
   a modified graph that directly incorporates spatial geometry into graph
   diffusion. See paper below for more. If spatial parameters aren't
   provided, this function reverts to a standard KNN graph.

        - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-
          Spatial Diffusion Geometry for Hyperspectral Image Clustering. 
          (In Review)

- If the NEigs is not included in the Hyperparameters structure, m = argmax( abs(diff(EigenVals)))
  eigenvalues are included in the Graph structure.

Inputs: X:                      Data matrix.
        Hyperparameters:        Structure with graph parameters with the 
                                required fields: 
                                    - Sigma>0: diffusion scale.
                                    - DiffusionNN: no. nearest neighbors.
                                    - NEigs: no. eigenvectors. (Optional)
        Dist:                   Matrix encoding pairwise distances between 
                                points in X. (Optional)

Output: Graph structure with the following fields:
            - Hyperparameters:  Parameters used to generate Graph structure.
            - W:                Weight matrix.
            - P:                Transition matrix for Markov diffusion
                                process.
            - StationaryDist:   Stationary distribution of P.
            - EigenVecs:        First m eigenvectors of P.
            - EigenVals:        First m eigenvalues of P.

© 2021 Sam L Polk, Tufts University. 
email: samuel.polk@tufts.edu
%}


% Extract Graph Parameters
sigma = Hyperparameters.Sigma;
NN = Hyperparameters.DiffusionNN;
n = length(X);

if nargin<3
    Dist=squareform(pdist(X)); %Pairwise distances between points in X
end


if isfield(Hyperparameters, 'SpatialParams') 

    % If spatial information is included in hyperparameters structure, we
    % incorporate that into the diffusion process.
    R = Hyperparameters.SpatialParams.SpatialRadius;
    M = Hyperparameters.SpatialParams.ImageSize(1);
    N = Hyperparameters.SpatialParams.ImageSize(2);
    
    % Build spatial nearest neighbors cells
    X_spatial=reshape(X,M,N,size(X,2));
    
    NN_Idx=cell(M,N);
    NN_Counts=zeros(M,N);
    NN_Dists=cell(M,N);
    map = zeros(M*N,2);
    for i=1:M
        for j=1:N
            NN_Idx{i,j}=FindNeighbors([i,j], R, M, N);
            NN_Counts(i,j)=length(NN_Idx{i,j});
            DistTemp=sqrt(sum((repmat(squeeze(X_spatial(i,j,:))',size(NN_Idx{i,j},2),1)-X(sub2ind([M,N], NN_Idx{i,j}(1,1:size(NN_Idx{i,j},2)), NN_Idx{i,j}(2,1:size(NN_Idx{i,j},2))),:)).^2,2))';
            [~,Idx]=sort(DistTemp,'ascend');
            NN_Idx{i,j}=NN_Idx{i,j}(:,Idx);
            NN_Idx{i,j}=sub2ind([M,N],NN_Idx{i,j}(1,:),NN_Idx{i,j}(2,:));% Convert back to subscript
            NN_Dists{i,j}=DistTemp(Idx);
            
            Temp = zeros(M,N);
            Temp(i,j) = 1;
            map(reshape(Temp,M*N,1) == 1,:) = [i,j];

        end
    end
        
    % Create weight matrix W
    W = zeros(n);
    P = zeros(n);
    D = zeros(n);
    for i = 1:n    
        
        i_row = map(i,1);
        i_col = map(i,2);
        
        if NN_Counts(i_row, i_col)>=NN
            D_sorted = NN_Dists{i_row,i_col}(1:NN);
            Idx = NN_Idx{i_row,i_col}(1:NN);
        else
            D_sorted = NN_Dists{i_row,i_col};
            Idx = NN_Idx{i_row,i_col};
        end
        W(i,Idx) = exp(-(D_sorted.^2)./(sigma^2));
        
        D(i,i) = sum(W(i,Idx));
        P(i,Idx) =  W(i,Idx)./D(i,i);
    end

else
    % Case where there is no spatial regularization. So, we build a standard KNN graph
        
    % Create weight matrix W, degree matrix D, and transition matrix P
    W = zeros(n);
    P = zeros(n);
    D = zeros(n);
    for i = 1:n    
        [D_sorted, sorting] = mink(Dist(i,:), NN+1); 
        W(i,sorting(2:end)) = exp(-(D_sorted(2:end).^2)./(sigma^2));
        D(i,i) = sum(W(i,:));
        P(i,sorting(2:end)) =  W(i,sorting(2:end))./D(i,i);
    end
end

% Calculate pi, stationary distribution
pi = diag(D)./sum(diag(D));

% Calculate the eigendecomposition of P
% try 
%     
%     [EigenVecs,EigenVals] = eigs(P, 11);
%     [EigenVals,idx] = sort(diag(abs(EigenVals)), 'descend');
%  
%     n_eigs = 11;
%     EigenVecs = real(EigenVecs(:,idx(1:n_eigs)));
%     EigenVals = EigenVals(1:n_eigs);
    
try 
 
    % If the number of eigenvalues is specified, choose that.
    if isfield(Hyperparameters, 'NEigs') 
        n_eigs = Hyperparameters.NEigs;
        [EigenVecs,EigenVals] = eigs(P, n_eigs);
        [EigenVals,idx] = sort(diag(abs(EigenVals)), 'descend');
           
    else
        % Otherwise, use the maximizer of the eigengap of P.

        [EigenVecs,EigenVals] = eigs(P, 20);
        [EigenVals,idx] = sort(diag(abs(EigenVals)), 'descend');
        [~,n_eigs] = max( abs(diff(EigenVals)));
        
        % Take care of fringe cases where we don't gather enough eigenvectors
        if n_eigs <5
            n_eigs = 5;
        end
    end
    EigenVecs = real(EigenVecs(:,idx(1:n_eigs)));
    EigenVals = EigenVals(1:n_eigs);


    % Set theoretical value for first eigenpair
    EigenVecs(:,1) = 1;
    EigenVals(1) = 1;

    % Store in graph structure
    Graph.Hyperparameters.Sigma = sigma;
    Graph.Hyperparameters.DiffusionNN = NN;
    Graph.EigenVecs = EigenVecs;
    Graph.EigenVals = EigenVals;
    Graph.StationaryDist = pi;
    Graph.P = P;
    Graph.W = W;
    
catch
    disp('EigenDecomposition of P failed')
    
    Graph.EigenVecs = NaN;
    Graph.EigenVals = [NaN, NaN];
    Graph.StationaryDist = NaN;
    Graph.P = NaN;
    Graph.W = NaN;
end


function Neighbors=FindNeighbors(Coordinates, Size, m, n)
%For a pixel coordinates in a matrix, find the size^2-1
%spatial nearest neighbors.

Neighbors=[];

for i=1:2*Size+1
    for j=1:2*Size+1
        if ( (Coordinates(1)-Size-1+i)>0 && (Coordinates(2)-Size-1+j)>0 && (Coordinates(1)-Size-1+i)<(m+1) && (Coordinates(2)-Size-1+j)<(n+1))
            Neighbors=[Neighbors,[Coordinates(1)-Size-1+i,Coordinates(2)-Size-1+j]'];
        end
    end
end


