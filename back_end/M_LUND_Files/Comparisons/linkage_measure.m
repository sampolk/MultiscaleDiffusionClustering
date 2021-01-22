function [Lin, Lbtw] = linkage_measure(X,C)

D = squareform(pdist(X)); % store distances beforehand

K = length(unique(C));
candidates_in = zeros(K,1);
candidates_btw = zeros(K,K);

for k=1:K
    candidates_in(k)  = max(D(C==k, C==k), [], 'all');
    
    cand_min = min(D(C==k,~(C==k)), [], 'all');
    
    if isempty(cand_min)
        cand_min = NaN;
    end
    candidates_btw(k) = cand_min;
end

Lin = max(candidates_in);
Lbtw = min(candidates_btw, [], 'all');
