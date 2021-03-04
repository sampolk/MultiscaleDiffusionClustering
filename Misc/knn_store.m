function [X,M,N,varargout] = knn_store(HSI, NN_max,save_on, noise_on)
%{
Objectives:     1. Standardize the points in a hyperspectral image (HSI).
                2. Compute nearest neighbor searches for each point (Optional).        
                3. Save to hard drive (Optional). 

Inputs:         HSI:        MxNxD HSI. 
                NN_max:     Number of nearest neighbors to search for.
                            (Optional)
                save_on:    Binary variable indicating whether to save. 
                            Results are saved if save_on == 1; (Optional)

                noise_on:   Binary variable indicating whether to add 
                            Gaussian noise. (Optional)

Outputs:        X:          (M*N)xD data matrix with standardized pixel spectra.
                M:          Number of row pixels.
                N:          Number of column pixels.
                
If NN_max is included in the input arguments, Dist_NN and Idx_NN are also
outputted:  
        Dist_NN:                nxM matrix where Dist_NN(i,:) encodes the M 
                                nearest neighbors of X(i,:), sorted in 
                                ascending order. M must be greater than NN.
                                (Optional)
        Idx_NN:                 nxM matrix where Idx_NN(i,:) encodes the 
                                indices of the M nearest neighbors of 
                                X(i,:) in X.  M must be greater than NN.
                                (Optional)

If save_on is not included, no save is made. 

If noise_on is not included, no noise is added to the HSI.

%}

if nargin <4
    noise_on = 1;
end

[M,N,~] = size(HSI); %dimensions of HSI
n = M*N; % number of pixels

% reshape and normalize data
X=reshape(HSI,size(HSI,1)*size(HSI,2),size(HSI,3));
if noise_on
    X = X + 10^(-6).*randn(size(X));
end
X=X./repmat(sqrt(sum(X.*X,1)),size(X,1),1); % Normalize HSI

if nargin>1
    % store distances 
    Dist_NN = zeros(n,NN_max);
    Idx_NN = zeros(n,NN_max);
    parfor i = 1:n
        [idx, dist_idx] = knnsearch(X, X(i,:), 'K', (NN_max+1));
        Dist_NN(i,:) = dist_idx(2:end);
        Idx_NN(i,:) = idx(2:end);
        disp(strcat('Nearest neighbors, ', num2str((1-i/n)*100, 3), '% complete.'))
    end
    varargout{1} = Dist_NN;
    varargout{2} = Idx_NN;
    
end

% save nearest neighbor searches and normalized X.
if nargin == 3
    if save_on == 1
        save('normalized_data_NNs.mat', 'X', 'M', 'N', 'Dist_NN', 'Idx_NN')
    end
end