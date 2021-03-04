function y = plot_abundance(X)

figure;
subplot(2,4,1)
imagesc(X(:,:,1))
axis equal off
title('EM1')

subplot(2,4,2)
imagesc(X(:,:,2))
axis equal off
title('EM2')

subplot(2,4,3)
imagesc(X(:,:,3))
axis equal off
title('EM3')

subplot(2,4,4)
imagesc(X(:,:,4))
axis equal off
title('EM4')

subplot(2,4,5)
imagesc(X(:,:,5))
axis equal off
title('EM5')

subplot(2,4,6)
imagesc(X(:,:,6))
axis equal off
title('EM6')

subplot(2,4,7)
imagesc(X(:,:,7))
axis equal off
title('EM7')
y = 1;

end