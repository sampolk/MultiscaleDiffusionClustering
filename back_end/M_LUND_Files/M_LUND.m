function Clusterings = M_LUND(X, Hyperparameters, G, p)
%{
 - This function produces a structure with multiscale clusterings produced
   with the M-LUND algorithm, presented in the following paper. 

        - Murphy, James M and Polk, Sam L., 2021. A Multiscale Environment 
          for Learning By Diffusion. (In Preparation)

   and analyzed further in the following paper:

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
                                Labels(:,t) reflects the LUND clustering at
                                time t.
            - K                 (T+2)x1 vector storing the number of 
                                clusters. K(t) reflects the number of
                                clusters in the Labels(:,t) clustering.
            - TotalVI           (T+2)x1 vector storing the Total VI of
                                clusterings. TotalVI(t) reflects the Total 
                                VI of the Labels(:,t) clustering.      
            - TimeSamples:      (T+2)x1 vector storing time steps during 
                                which LUND was evaluated. TimeSamples(t) 
                                reflects the time step at which LUND was
                                evaluted to generate Labels(:,t).
            - Density:          Kernel Density Estimator used in LUND
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
        [Ct(:,i),Kt(i), Dt(:,i)] = LearningbyUnsupervisedNonlinearDiffusion(X, timesamples(i), G, p);
    end

    % ============================ VI analysis ============================

    J = find(and(Kt<n/2, Kt>1)); % Time samples during which a nontrivial clusering is extracted
    n_J = length(J);

    if n_J > 0 % There is a nontrivial clustering of X.
        [~,VI_tot,t] = totalVI_minimization(Ct, Kt);
        TotalVI.Vector = VI_tot;
        TotalVI.Minimizer_Idx = t;
    else
        TotalVI.Vector = zeros(T+2,1);
        TotalVI.Vector(:) = NaN; 
        TotalVI.Minimizer_Idx = NaN;
    end

    % Store results in structure "Clusterings"
    Clusterings.Graph = G;
    Clusterings.Hyperparameters = Hyperparameters;
    Clusterings.Labels = Ct;
    Clusterings.K = Kt;
    Clusterings.TotalVI = TotalVI;
    Clusterings.TimeSamples = timesamples;
    Clusterings.Density = p;
    Clusterings.Dt = Dt;

else
    Clusterings = NaN;
end