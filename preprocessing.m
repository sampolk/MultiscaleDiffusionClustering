%% Load data and process

load('r1_reg.mat')

HSI = r1_reg;
NN_max = 3500;
save_on = 0;
noise_on = 1;

[X,M,N,Dist_NN, Idx_NN] = knn_store(HSI, NN_max,save_on, noise_on);
D1 = Dist_NN(:,1:1000);
D2 = Dist_NN(:,1001:2000);
D3 = Dist_NN(:,2001:3000);
D4 = Dist_NN(:,3001:3500);
Idx1 = Idx_NN(:,1:1000);
Idx2 = Idx_NN(:,1001:2000);
Idx3 = Idx_NN(:,2001:3000);
Idx4 = Idx_NN(:,3001:3500);
clear 'Dist_NN' 'Idx_NN'
save('normalized_data_NNs.mat')