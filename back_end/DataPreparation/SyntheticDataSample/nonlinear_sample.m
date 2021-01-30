function [X,Y] = nonlinear_sample(rad_in, rad_mid, rad_out, n_in, n_mid, n_out, max_spl)

flag = 0;
% Finds a suitable sample
for k = 1:max_spl
    if flag == 0

        C1 = (rad_in.*(0.1*randn(n_in,1)+1)).*exp(1i.*rand(n_in,1)*2*pi);
        C2 = (rad_mid.*(0.04*randn(n_mid,1)+1)).*exp(1i.*rand(n_mid,1)*2*pi);
        C3 = (rad_out.*(0.02*randn(n_out,1)+1)).*exp(1i.*rand(n_out,1)*2*pi);

        X = [C1; C2; C3];
        X = [real(X), imag(X)];
        Y = [ones(n_in,1); 2.*ones(n_mid,1); 3.*ones(n_out,1)];

        load('nonlinear-HP.mat')
        G = extract_graph(X, Hyperparameters);

        if G.EigenVals(2)<1 

            Clusterings = M_LUND(X, Hyperparameters);

            n_K_nt = length(unique(Clusterings.K(find(and(Clusterings.K<length(X)/2, Clusterings.K>1)))));
            
            if length(intersect(unique(Clusterings.K), [2,3])) == 2 && n_K_nt==2
                
                flag = 1;

            end

        end
    end
end