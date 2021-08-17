function [C,Ks, t_vals, Alpha, Beta] = HierearchicalSpectralClustering(G, T_max)

lambda = sort(G.EigenVals, 'descend');
[~,max_size] = computer;

% First, check to see if RAM requirements are too high. 
if T_max >= floor(max_size/10^8) % 1e-7 times the maximum size allocated in MATLAB

    % In this case, we do exponential sampling of the time domain. 
    base = T_max^(1./(floor(max_size/10^8)-1));
    TimeSamples = base.^(0:floor(max_size/10^8)-1);
    
    Deltat = zeros(length(TimeSamples),1);
    Kt = zeros(length(TimeSamples),1);
    for t = 1:length(TimeSamples)
        [Deltat(t), Kt(t)] = max(lambda(1:end-1).^TimeSamples(t) -lambda(2:end).^TimeSamples(t));
    end
    
else
    % In this case, we do linear sampling of the time domain.
    Deltat = zeros(T_max,1);
    Kt = zeros(T_max,1);
    for t = 1:T_max
        [Deltat(t), Kt(t)] = max(lambda(1:end-1).^t -lambda(2:end).^t);
    end
end

% Use MATLAB findpeaks built-in to find the t maxima. 
[Delta, Ts] = findpeaks(Deltat);

[d,t]= max(Deltat);
if isempty(intersect(Ts,t))
    Delta = [Delta; d]; 
    Ts = [Ts; t];
end
Ks = Kt(Ts);
Ts = [0; Ts]; % Add t0 = 0.

% Initialize Cluster analysis
M = length(Delta); % No. unique clusterings
C = zeros(length(G.EigenVecs), M);
Alpha = zeros(M,1);
t_vals = zeros(M,2);
Beta = zeros(M,1);

% Cluster M clusterings based on HSC scheme 
parfor l = 1:M
    C(:,l) = SpectralClustering(Ks(l), G);
    t_vals(l,:) = Ts(l:l+1);
    Alpha(l) = (Ts(l+1) - Ts(l))/T_max;
    Beta(l) = Delta(l);
end
