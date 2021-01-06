function [X,Y] = bottleneck_sample(n_g, n_u, std_g,  width_u, height_u, separation)


% ============================ Left Bottleneck ============================
% Bottleneck 1 Gaussians
gaussian_up = std_g.*zscore(randn(n_g, 2)) + [0,height_u] + [0, std_g]; % Gaussian on top
gaussian_low = std_g.*zscore(randn(n_g, 2)) - [0, std_g]; % Gaussian on bottom

% Delete outliers in Gaussians

gaussian_up_norm = gaussian_up- mean(gaussian_up);
gaussian_up(or(abs(gaussian_up_norm(:,1))>2*std_g, abs(gaussian_up_norm(:,2))>2*std_g),:) = [];

gaussian_low_norm = gaussian_low- mean(gaussian_low);
gaussian_low(or(abs(gaussian_low_norm(:,1))>2*std_g, abs(gaussian_low_norm(:,2))>2*std_g),:) = [];

% Uniform distribution
uniform(:,1) = width_u.*rand(n_u,1)-0.5;
uniform(:,2) = height_u.*rand(n_u,1);

bottleneck1 = [gaussian_up; uniform; gaussian_low];

% =========================== Right Bottleneck ============================

% Bottleneck 1 Gaussians
gaussian_up = std_g.*zscore(randn(n_g, 2)) + [0,height_u] + [0, std_g]; % Gaussian on top
gaussian_low = std_g.*zscore(randn(n_g, 2)) - [0, std_g]; % Gaussian on bottom

% Delete outliers in Gaussians

gaussian_up_norm = gaussian_up- mean(gaussian_up);
gaussian_up(or(abs(gaussian_up_norm(:,1))>2*std_g, abs(gaussian_up_norm(:,2))>2*std_g),:) = [];

gaussian_low_norm = gaussian_low- mean(gaussian_low);
gaussian_low(or(abs(gaussian_low_norm(:,1))>2*std_g, abs(gaussian_low_norm(:,2))>2*std_g),:) = [];

% Uniform distribution
uniform(:,1) = width_u.*rand(n_u,1)-0.5;
uniform(:,2) = height_u.*rand(n_u,1);

bottleneck2 = [gaussian_up; uniform; gaussian_low] + [separation,0];

% ========================== Separated Gaussian ===========================

gaussian_sep = std_g.*zscore(randn(n_g, 2)) + [0.7*separation,height_u + std_g] + [0, 4*std_g]; % Separated Gaussian

gaussian_sep_norm = gaussian_sep- mean(gaussian_sep);
gaussian_sep(or(abs(gaussian_sep_norm(:,1))>2*std_g, abs(gaussian_sep_norm(:,2))>2*std_g),:) = [];

% ============================= Collect Data  =============================

X = [bottleneck1; gaussian_sep; bottleneck2];
Y = [ones(length(bottleneck1),1); 2*ones(length(gaussian_sep),1); 3*ones(length(bottleneck2),1)];


