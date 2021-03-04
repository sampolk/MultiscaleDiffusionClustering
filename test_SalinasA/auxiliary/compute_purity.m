function purity = compute_purity(X,h,w,k)
% Input: X - 3D Abundance map
%        h - rows
%        w - columns
%        k - number of endmembers
% Output: Highest abundance
purity = reshape(max(reshape(X, h*w,k)')', h,w);
figure;
imagesc(purity)
axis equal off

purity = reshape(purity, h*w, 1);

end