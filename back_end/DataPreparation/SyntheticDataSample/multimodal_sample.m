function [X,Y] = multimodal_sample(nGaussians, nUniform)

if nargin == 0
    nGaussians = 100; 
    nUniform = 50;
end

sigma1 = 0.64;
sigma2 = 0.11;
sigma3 = 0.11;

% 1st and 4th quadrant
X1 = [(1+0.01*randn((6*nGaussians),1)).*exp(1i*( sigma1*randn((6*nGaussians),1))); (1+0.01*randn(floor(4*nUniform),1)).*exp(1i*( pi*rand(floor(4*nUniform),1)-pi/2))];

% 2nd quadrant
X2 = [(1+0.01*randn((3*nGaussians),1)).*exp(1i*( sigma2*randn((3*nGaussians),1)+(5/8+1/32)*pi)); (1+0.01*randn(floor(nUniform),1)).*exp(1i*( pi/4*rand(floor(nUniform),1)+pi/2))];
X3 = [(1+0.01*randn((3*nGaussians),1)).*exp(1i*( sigma2*randn((3*nGaussians),1)+(7/8-1/32)*pi)); (1+0.01*randn(floor(nUniform),1)).*exp(1i*( pi/4*rand(floor(nUniform),1)+3*pi/4))];

% 3rd quadrant
X4 = [(1+0.01*randn((3*nGaussians),1)).*exp(1i*( sigma3*randn((3*nGaussians),1)+(9/8+1/32)*pi)); (1+0.01*randn(floor(nUniform),1)).*exp(1i*( pi/4*rand(floor(nUniform),1)+pi))];
X5 = [(1+0.01*randn((3*nGaussians),1)).*exp(1i*( sigma3*randn((3*nGaussians),1)+ (11/8-1/32)*pi)); (1+0.01*randn(floor(nUniform),1)).*exp(1i*( pi/4*rand(floor(nUniform),1)+5*pi/4))];

X = [X1; X2; X3; X4; X5];
X = [real(X), imag(X)];

Y = [ones(size(X1)); 2*ones(size(X2)); 3*ones(size(X3));  4*ones(size(X4)); 5*ones(size(X5))];
