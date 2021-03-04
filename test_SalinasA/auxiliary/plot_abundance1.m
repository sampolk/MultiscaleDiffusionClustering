function y = plot_abundance1(X)

figure;
subplot(2,3,1)
imagesc(X(:,:,1))
axis equal off
title('EM1')

subplot(2,3,2)
imagesc(X(:,:,2))
axis equal off
title('EM2')

subplot(2,3,3)
imagesc(X(:,:,3))
axis equal off
title('EM3')

subplot(2,3,4)
imagesc(X(:,:,4))
axis equal off
title('EM4')

subplot(2,3,5)
imagesc(X(:,:,5))
axis equal off
title('EM5')

subplot(2,3,6)
imagesc(X(:,:,6))
axis equal off
title('EM6')
y = 1;

end