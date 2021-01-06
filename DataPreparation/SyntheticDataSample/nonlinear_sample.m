function [X,Y] = nonlinear_sample(rad_in, rad_mid, rad_out, n_in, n_mid, n_out)

C1 = (rad_in.*(0.1*randn(n_in,1)+1)).*exp(1i.*rand(n_in,1)*2*pi);
C2 = (rad_mid.*(0.04*randn(n_mid,1)+1)).*exp(1i.*rand(n_mid,1)*2*pi);
C3 = (rad_out.*(0.02*randn(n_out,1)+1)).*exp(1i.*rand(n_out,1)*2*pi);

X = [C1; C2; C3];
X = [real(X), imag(X)];
Y = [ones(n_in,1); 2.*ones(n_mid,1); 3.*ones(n_out,1)];