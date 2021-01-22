function [delta, lambda, kappa] = StochasticComplement(P,C)

K = length(unique(C));
n = length(P);
S = eye(n); % initialize with the identity.

for k = 1:K
    within_idces = find(C == k);
    if length(within_idces)>1 % Not singleton
        outside_idces = setxor(1:n, within_idces);

        Pkk = P(within_idces, within_idces);
        Pk = P(outside_idces, outside_idces);
        PkStar = P(within_idces, outside_idces);
        PStark = P(outside_idces, within_idces);
        
        S(within_idces, within_idces) = Pkk + PkStar*((eye(length(Pk)) - Pk)\PStark);        
    end
end

[Z,D] = eig(S);
EigenVals = sort(abs(diag(D)), 'descend');
lambda = EigenVals(K+1);
delta = norm(P-S, inf);
kappa = cond(Z, inf);