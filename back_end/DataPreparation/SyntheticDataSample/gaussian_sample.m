function [X,Y] = gaussian_sample(std_in, std_out, n_in, n_out, max_spl)


flag = 0;
% Finds a suitable sample
for k = 1:max_spl
    if flag == 0
        % Sample data
        X_in_pos = randn(n_in,3);
        X_in_neg = randn(n_in,3);
        X_out_pos = randn(n_out, 3);
        X_out_neg = randn(n_out, 3);

        % Standardize and implement mean and standard deviations as needed.
        X_in_pos = std_in.*zscore(X_in_pos) + [3*std_in, 0,0];
        X_in_neg = std_in.*zscore(X_in_neg) - [3*std_in, 0,0];
        X_out_pos = std_out.*zscore(X_out_pos) + [6*std_in, 0,0] + [2*std_out,0,0];
        X_out_neg = std_out.*zscore(X_out_neg) - [6*std_in, 0,0] - [2*std_out,0,0];

        Y = [ones(n_out,1); 2*ones(n_in,1); 3*ones(n_in,1); 4*ones(n_out,1)];

        X = [ X_out_neg; X_in_neg; X_in_pos; X_out_pos];

        load('gaussians-HP.mat')
        G = extract_graph(X, Hyperparameters);

        if G.EigenVals(2)<1 

            Clusterings = M_LUND(X, Hyperparameters);
            

            if length(intersect(unique(Clusterings.K), [2,4])) == 2  
                
                if sum(Clusterings.K ==4)>1 && sum(Clusterings.K ==2)>1 && Clusterings.K(1)>=4
                    flag = 1;
                end
            end

        end
    end
end
