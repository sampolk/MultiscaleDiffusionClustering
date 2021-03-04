%% Pre-process the data
clear
[X,Y] = extract_salinasA();
data_name = 'SalinasA';
load('salinasA-HP.mat')

D = squareform(pdist(X));
dist = D(D>0);

G = extract_graph(X, Hyperparameters, D);


%% ATGP + LUND, K = 7
k=length(unique(Y));
X_col = X';
T = 11:15;
[acc_atgp_7,labels_atgp_7,t_atgp_7] = atgp_lund(X_col,k,83,86,G,T,Y,0);


%% ATGP + LUND, K = 6
k=length(unique(Y))-1;
X_col = X';
T = 5:12;
[acc_atgp_6,labels_atgp_6,t_atgp_6] = atgp_lund(X_col,k,83,86,G,T,Y,1);


%% N-Finder + LUND, K = 7
k=length(unique(Y));
X_col = X';
T = 6:12;
[acc_nfindr_7,labels_nfindr_7,t_nfindr_7] = nfindr_lund(X_col,k,83,86,G,T,Y,0);


%% N-Finder + LUND, K = 6
k=length(unique(Y))-1;
X_col = X';
T = 5:11;
[acc_nfindr_6,labels_nfindr_6,t_nfindr_6] = nfindr_lund(X_col,k,83,86,G,T,Y,1);


%% N-Finder-I + LUND, K = 7
k=length(unique(Y));
X_col = X';
T = 5:11;
[acc_nfindri_7,labels_nfindri_7,t_nfindri_7] = nfindri_lund(X_col,k,83,86,G,T,Y,0);


%% N-Finder-I + LUND, K = 6
k=length(unique(Y))-1;
X_col = X';
T = 6:12;
[acc_nfindri_6,labels_nfindri_6,t_nfindri_6] = nfindri_lund(X_col,k,83,86,G,T,Y,1);


%% VCA + LUND, K = 7
k=length(unique(Y));
X_col = X';
T = 5:11;
[acc_vca_7,labels_vca_7,t_vca_7] = vca_lund(X_col,k,83,86,G,T,Y,0);


%% VCA + LUND, K = 6
k=length(unique(Y))-1;
X_col = X';
T = 5:2:13;
[acc_vca_6,labels_vca_6,t_vca_6] = vca_lund(X_col,k,83,86,G,T,Y,1);


%% VCA-I + LUND, K = 7
k=length(unique(Y));
X_col = X';
T = 5:11;
[acc_vcai_7,labels_vcai_7,t_vcai_7] = vcai_lund(X_col,k,83,86,G,T,Y,0);


%% VCA-I + LUND, K = 6
k=length(unique(Y))-1;
X_col = X';
T = 5:11;
[acc_vcai_6,labels_vcai_6,t_vcai_6] = vcai_lund(X_col,k,83,86,G,T,Y,1);


%% KDE + LUND, K = 7
p = KDE(X, Hyperparameters, D);
figure
imagesc(reshape(p, 83,86))
title('p')
axis equal off

T = 12;
[acc_kde_7,labels_kde_7,~,t_kde_7] = LUND_labels(X, Y, p, G, T, 83, 86, 7, 0);


%% KDE + LUND, K = 6
p = KDE(X, Hyperparameters, D);
figure
imagesc(reshape(p, 83,86))
title('p')
axis equal off

T = 5;
[acc_kde_6,labels_kde_6,~,t_kde_6] = LUND_labels(X, Y, p, G, T, 83, 86, k, 1);


