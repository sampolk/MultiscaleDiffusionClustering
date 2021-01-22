function [V, NMI] = VI(U,V)
%{
Calculates the variation of information and normalized mutual information
between clusterings U and V of a dataset X. 

Inputs:     U:  nx1 vector encoding first clustering of X.
            V:  nx1 vector encoding second clustering of X.

Outputs:    V:  The variation of information between U and V. 
            U:  The normalized mutual information between U and V. 

%}
% See if dimensions are right
assert(numel(U) == numel(V));

% initialize
n = numel(U);
V = reshape(U,1,n);
V = reshape(V,1,n);

% ensures clusterings start at K=1. 
l = min(min(U),min(V));
U = U-l+1;
V = V-l+1;
K = max(max(U),max(V));

% Calculate the joint and  densities. 
idx = 1:n;
Mu = sparse(idx,U,1,n,K,n);
Mv = sparse(idx,V,1,n,K,n);
Puv = nonzeros(Mu'*Mv/n);       % joint probability distribution of U and V
Huv = -dot(Puv,log2(Puv));      % joint entropy of U and V
Pu = nonzeros(mean(Mu,1));      % probability distribution of clustering U.
Pv = nonzeros(mean(Mv,1));      % probability distribution of clustering V.

Hu = -dot(Pu,log2(Pu));         % entropy of Pu
Hv = -dot(Pv,log2(Pv));         % entropy of Pv
MI = Hu + Hv - Huv;             % mutual information between Pu & Pv

% normalized mutual information between Pu & Pv
NMI = sqrt((MI/Hu)*(MI/Hv));    
NMI = MI/(Hu+Hv)*2;
NMI = max(0,NMI);

V = Hu + Hv - 2*MI; 
