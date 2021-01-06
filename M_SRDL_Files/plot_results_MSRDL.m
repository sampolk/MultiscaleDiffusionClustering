function fig = plot_results_MSRDL(Clusterings_SRDL, Clusterings_LUND)

n = length(Clusterings_SRDL.Labels);
M = Clusterings_SRDL.Hyperparameters.SpatialParams.ImageSize(1);
N = Clusterings_SRDL.Hyperparameters.SpatialParams.ImageSize(2);

% -------------------------- Find K for M-SRDL ----------------------------
J = find(and(Clusterings_SRDL.K>1, Clusterings_SRDL.K<n/2));
K_nt_srdl = unique(Clusterings_SRDL.K(J));    

K_nt_srdl = K_nt_srdl(end:-1:1);

if Clusterings_SRDL.K(1)>n/2 || Clusterings_SRDL.K(1) <1
    K_nt_srdl = [Clusterings_SRDL.K(1); K_nt_srdl];
end
if Clusterings_SRDL.K(end)>n/2 || Clusterings_SRDL.K(end) <1
    K_nt_srdl = [K_nt_srdl; Clusterings_SRDL.K(end)];
end

% -------------------------- Find K for M-LUND ----------------------------

J = find(and(Clusterings_LUND.K>1, Clusterings_LUND.K<n/2));
K_nt_lund = unique(Clusterings_LUND.K(J));    

% fig = figure;
K_nt_lund = K_nt_lund(end:-1:1);
n_plots = length(K_nt_lund);

if Clusterings_LUND.K(1)>n/2 || Clusterings_LUND.K(1) <1
    K_nt_lund = [Clusterings_LUND.K(1); K_nt_lund];
end
if Clusterings_LUND.K(end)>n/2 || Clusterings_LUND.K(end) <1
    K_nt_lund = [K_nt_lund; Clusterings_LUND.K(end)];
end

fig = figure;
n_cols = max([length(K_nt_srdl), length(K_nt_lund)]);


% -------------------------- Plot M-SRDL Results --------------------------
row = 1;
for col = 1:n_cols
    
    subplot( 2, n_cols, (row-1)*n_cols + col)
    
    if col-1 < n_cols
        t = find(Clusterings_SRDL.K == K_nt_srdl(col), 1, 'first');
    else
        t = length(Clusterings_SRDL.K);
    end
    imagesc(reshape(Clusterings_SRDL.Labels(:,t), M,N))
    title(strcat('SRDL Clustering, K=', num2str(K_nt_srdl(col))), 'interpreter', 'latex')
    xticks([])
    yticks([])
    pbaspect([1,1,1])
    set(gca,'FontSize', 20, 'FontName', 'Times')
end

% -------------------------- Plot M-LUND Results --------------------------
row = 2;
for col = 1:n_cols
    
    subplot( 2, n_cols, (row-1)*n_cols + col)
    
    if col < n_cols
        t = find(Clusterings_LUND.K == K_nt_lund(col), 1, 'first');
    else
        t = length(Clusterings_LUND.K);
    end
    imagesc(reshape(Clusterings_LUND.Labels(:,t), M,N))
    title(strcat('LUND Clustering, K=', num2str(K_nt_lund(col))), 'interpreter', 'latex')
    xticks([])
    yticks([])
    pbaspect([1,1,1])
    set(gca,'FontSize', 20, 'FontName', 'Times')
end