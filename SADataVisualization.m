%{
This script loads the Salinas A dataset and visualizes its ground truth
labels and a random sample of its pixels' spectra. 

This figure appears as in the following papers:

     - Murphy, James M and Polk, Sam L. "A Multiscale Environment for 
       Learning By Diffusion." In Preparation (2021).
     - Polk, Sam L. and Murphy James M. "Multiscale Spectral-Spatial 
       Diffusion Geometry for Hyperspectral Image Clustering." In Review 
       (2021).

(c) Sam L. Polk: samuel.polk@tufts.edu
%}

%% Extract data
[X,Y] = extract_salinasA();
M = 83;
N = 86;
GT = reshape(Y, M,N);

%% Figure 6

figure

colors = {'k', '#4658F8', '#2896EB', '#13BEB8', '#80CA57', '#FCBB3D', '#F8FA13'};
rgb_key = [[0,0,0]; [0.275,0.345,0.973]; [15.7,58.8, 92.2]./100; [7.5, 74.5, 72.2]./100; [50.2,79.2,34.1]./100; [98.8, 73.3, 23.9]./100; [97.3, 98,7.5]./100];

c_data = zeros(M,N, 3);
for i = 1:M
    for j = 1:N
        c_data(i,j,:) = rgb_key(GT(i,j),:);
    end
end

subplot(1,2,1)
image(c_data)
title('Salinas A Ground Truth', 'interpreter', 'latex')
xticks([])
yticks([])
pbaspect([1,1,1])
set(gca,'FontName', 'Times', 'FontSize', 20)

subplot(1,2,2)
sample = randsample(length(X),150);
hold on
for i = 1:length(sample)
    c = colors(Y(sample(i)));
    plot(1:224, X(sample(i),:),'Color', c{1})
end
box on
title('Randomly Selected Pixel Spectra, Colored by Class', 'interpreter', 'latex')
xlabel('Spectral Band Number')
ylabel('Recorded Reflectance')

pbaspect([1,1,1])
set(gca,'FontName', 'Times', 'FontSize', 20)
