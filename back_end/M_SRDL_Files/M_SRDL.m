function Clusterings = M_SRDL(X, Hyperparameters, G, p)
%{
 - This function produces a structure with multiscale clusterings produced
   with the M-SRDL algorithm, presented in the following paper. 

        - Polk, Sam L. and Murphy James M., 2021. Multiscale Spectral-
          Spatial Diffusion Geometry for Hyperspectral Image Clustering. 
          (In Review)

Inputs: X:                      Data matrix.
        Hyperparameters:        Structure with graph parameters with the 
                                required fields: 
                                    - Sigma>0: diffusion scale.
                                    - DiffusionNN: no. nearest neighbors in
                                                   KNN graph.
                                    - NEigs: no.   eigenvectors. (Optional)
                                    - Sigma0>0:    KDE Bandwidth.
                                    - DesnityNN:   no. nearest neighbors  
                                                   used in KDE.
        G:                      Graph structure computed using  
                                'extract_graph.m' (Optional).
        p:                      Kernel Density Estimator (Optional).

Output: Clusterings structure with the following fields:

            - Graph:            Graph structure computed with
                                'extract_graph.m
            - Hyperparameters:  Parameters used to generate Graph structure.
            - Labels:           n x (T+2) matrix storing clusterings of X.
                                Labels(:,t) reflects the SRDL clustering at
                                time t.
            - K                 (T+2)x1 vector storing the number of 
                                clusters. K(t) reflects the number of
                                clusters in the Labels(:,t) clustering.
            - TotalVI           (T+2)x1 vector storing the Total VI of
                                clusterings. TotalVI(t) reflects the Total 
                                VI of the Labels(:,t) clustering.      
            - TimeSamples:      (T+2)x1 vector storing time steps during 
                                which SRDL was evaluated. TimeSamples(t) 
                                reflects the time step at which SRDL was
                                evaluted to generate Labels(:,t).
            - Density:          Kernel Density Estimator used in SRDL
                                clustering. nx1 vector.
            - Dt:               n x (T+2) matrix storing \mathcal{D}_t(x). 
                                Dt(:,t) reflects the values taken by
                                \mathcal{D}_t(x), where t = TimeSamples(t).

© 2021 Sam L Polk, Tufts University. 
email: samuel.polk@tufts.edu
%}

if nargin == 2
    G = extract_graph(X, Hyperparameters);
    p = KDE(X,Hyperparameters);
elseif nargin == 3
    p = KDE(X,Hyperparameters);
end

n = length(X);

T = full(ceil(log( log(Hyperparameters.Tau*min(G.StationaryDist)/2)/log(G.EigenVals(2)))/log(Hyperparameters.Beta)));

if isreal(T)
    
    % ========================== Cluster analysis =========================

    % Extract Time Steps
    timesamples = [0, Hyperparameters.Beta.^(0:T)];

    % Initialize
    Ct = zeros(n,T+2);
    Kt = zeros(T+2,1);
    Dt = zeros(n,T+2);
    parfor i = 1:T+2
        [Ct(:,i),Kt(i), Dt(:,i)] =SRDL(X, timesamples(i), Hyperparameters, G, p);
    end

    % ============================ VI analysis ============================

    J = find(and(Kt<n/2, Kt>1)); % Time samples during which a nontrivial clusering is extracted
    n_J = length(J);
    V = zeros(n_J); % To be the VI(Cs, Ct) matrix for nontrivial clusterings

    if n_J > 0 % There is a nontrivial clustering of X.
        notJ = setxor((1:T+2), J);
        for i = 1:n_J
            parfor j = 1:n_J
                V(i,j) = VI(Ct(:,J(i)), Ct(:,J(j)));
            end
        end
        % Calculate Total VI
        VI_tot = zeros(T+2,1);
        VI_tot(J) = sum(V);
        VI_tot(notJ) = NaN; % Assign NaN Total VI for trivial clusterings.
    else
        VI_tot = NaN;
    end

    % Store results in structure "Clusterings"
    Clusterings.Graph = G;
    Clusterings.Hyperparameters = Hyperparameters;
    Clusterings.Labels = Ct;
    Clusterings.K = Kt;
    Clusterings.TotalVI = VI_tot;
    Clusterings.TimeSamples = timesamples;
    Clusterings.Density = p;
    Clusterings.Dt = Dt;

else
    Clusterings = NaN;
end