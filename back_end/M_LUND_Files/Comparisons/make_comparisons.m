function [C_MMS, C_HSC, C_SLC, C_SLLUND] = make_comparisons(X, C_MLUND, mms_in, plot_on)

M = 83;
N = 86;

G = C_MLUND.Graph;
TimeSamples = C_MLUND.TimeSamples;
disp('M-LUND run complete')

% MMS Clustering
if mms_in

    [Stb, N1, ~, C1] = stability(G.W,TimeSamples,'full');

    C_MMS.Stability = Stb;
    C_MMS.K = N1;
    C_MMS.Labels = C1;
    C_MMS.TimeSamples = TimeSamples;
    C_MMS.Graph = G;

    disp('MMS Clustering run complete.')
else
    C_MMS = NaN;
end

% HSC
[C, K, t_vals, Alpha, Beta] = HierearchicalSpectralClustering(G, max(TimeSamples));

C_HSC.Graph = G;
C_HSC.Labels = C;
C_HSC.K = K;
C_HSC.Alpha = Alpha;
C_HSC.Beta = Beta;
C_HSC.TimeSamples = t_vals;
disp('HSC run complete.')

% SLC
Z_SLC = linkage(X,'single');
C_SLC  = ones(length(X),10);
for K = 1:5
    C_SLC(:,K) = cluster(Z_SLC,'MaxClust', K);
end
disp('SLC Clustering run complete.')

% SL-LUND
[~, C_SLLUND] = SL_LUND(X, C_MLUND);

if plot_on
    
    figure 
    
    i5 = find(C_MLUND.K == 5, 1,'first');
    i2 = find(C_MLUND.K == 2, 1,'first');
    
    % M-LUND plots
    subplot(2,4,1)
    imagesc(reshape(C_MLUND.Labels(:,i5),M,N))   
    title(strcat('M-LUND Clustering of Salinas A, K=5, log2(t)=',num2str(i5-2)))
    pbaspect([1,1,1])
    set(gca,'FontName', 'Times', 'FontSize', 12)
    xticks([])
    yticks([])    
    
    subplot(2,4,5)
    imagesc(reshape(C_MLUND.Labels(:,i2),M,N))    
    title(strcat('M-LUND Clustering of Salinas A, K=2, log2(t)=',num2str(i2-2)))
    pbaspect([1,1,1])
    set(gca,'FontName', 'Times', 'FontSize', 12)
    xticks([])
    yticks([])    

    i5 = find(C_HSC.K == 5, 1,'first');
    i2 = find(C_HSC.K == 2, 1,'first');
    
    % HSC plots
    subplot(2,4,2)
    imagesc(reshape(C_HSC.Labels(:,i5),M,N))
    title(strcat('HSC Clustering of Salinas A, K=5, log2(t)=',num2str(log2(t_vals(i5,1)))))
    pbaspect([1,1,1])
    set(gca,'FontName', 'Times', 'FontSize', 12)
    xticks([])
    yticks([])    
    
    subplot(2,4,6)
    imagesc(reshape(C_HSC.Labels(:,i2),M,N))    
    title(strcat('HSC Clustering of Salinas A, K=2, log2(t)=',num2str(log2(t_vals(i2,1)))))
    pbaspect([1,1,1])
    set(gca,'FontName', 'Times', 'FontSize', 12)
    xticks([])
    yticks([])        
    
    % SLC plots
    subplot(2,4,3)
    imagesc(reshape(C_SLC(:,5),M,N))    
    title(strcat('SLC Clustering of Salinas A, K=5'))

    pbaspect([1,1,1])
    set(gca,'FontName', 'Times', 'FontSize', 12)
    xticks([])
    yticks([])    
    
    subplot(2,4,7)
    imagesc(reshape(C_SLC(:,2),M,N))   
    title(strcat('SLC Clustering of Salinas A, K=2'))
    pbaspect([1,1,1])
    set(gca,'FontName', 'Times', 'FontSize', 12)
    xticks([])
    yticks([])   
    disp('M-LUND vs. HSC vs. SLC plot complete.')

    % SLC plots
    subplot(2,4,4)
    imagesc(reshape(C_SLLUND(:,1),M,N))    
    title(strcat('SL-LUND Clustering of Salinas A, K=5'))

    pbaspect([1,1,1])
    set(gca,'FontName', 'Times', 'FontSize', 12)
    xticks([])
    yticks([])    
    
    subplot(2,4,8)
    imagesc(reshape(C_SLLUND(:,4),M,N))   
    title(strcat('SL-LUND Clustering of Salinas A, K=2'))
    pbaspect([1,1,1])
    set(gca,'FontName', 'Times', 'FontSize', 12)
    xticks([])
    yticks([])   
    disp('M-LUND vs. HSC vs. SLC plot complete.')
    
    if mms_in
        
        Ks = unique(C_MMS.K(and(C_MMS.K<length(C_MMS.Labels)/2, C_MMS.K>1)));
        n_K = length(Ks);
        
        figure
        for k = 1:n_K
            K = Ks(k);
            i = find(C_MMS.K == K, 1, 'first');
            
            subplot(1,n_K,k)
            title(strcat('M-LUND Clustering of Salinas A, K=', num2str(K), ', log2(t)=',num2str(i-2)))
            imagesc(reshape(C_MLUND.Labels(:,i),M,N))    
            pbaspect([1,1,1])
            set(gca,'FontName', 'Times', 'FontSize', 12)
            xticks([])
            yticks([])   
    
        end
        disp('MMS plot complete.')

    end
end        

