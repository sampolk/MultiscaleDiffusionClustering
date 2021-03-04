function [acc,labels,CM,t_best, OA] = SRDL_labels(X, Y, p, Hyperparameters, G, T, h, w, k, flag)

t = [0,2.^(0:T)];
C = zeros(h*w,T+2);
OA = zeros(T+2,1);
Hyperparameters.K_known = k;

if flag == 1
    mask = Y;
    for i = 1:length(mask)
        if mask(i,1) ~= 1
            mask(i,1) = 0;
        end
    end
else
    mask = zeros(length(Y),1);
end





for i = 1:length(t)
    [C(:,i),~,~] = SRDL(X, t(i), Hyperparameters, G, p);
end

for i = 1:length(t)
    C_temp = C(:,i);
    C_t = C_temp(mask==0);
    Y_t = Y(mask==0);
    [OA(i), ~, ~] = calculateAccuracy(C_t, Y_t);
end

[~,I] = max(OA);
t_best = t(I);
[Cb,~,~] = SRDL(X, t_best, Hyperparameters, G, p);
[acc, labels, CM] = calculateAccuracy(Cb(mask==0), Y(mask==0));
labels1 = ones(length(Y),1);
labels1(mask==0) = labels;


figure;
subplot(1,2,1)
imagesc(reshape(Y,83,86))
title('GT')
axis equal off
set(gca,'FontSize', 20, 'FontName', 'Times')   
subplot(1,2,2)
%imagesc(reshape(C(:,I),83,86))
imagesc(reshape(labels1,83,86))
title('Endmembers+LUND')
axis equal off
set(gca,'FontSize', 20, 'FontName', 'Times') 