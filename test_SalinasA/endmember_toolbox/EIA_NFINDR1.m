% // ====================================================================
% // This file is part of the Endmember Induction Algorithms Toolbox for MATLAB 
% // Copyright (C) Grupo de Inteligencia Computacional, Universidad del 
% // País Vasco (UPV/EHU), Spain, released under the terms of the GNU 
% // General Public License.
% //
% // Endmember Induction Algorithms Toolbox is free software: you can redistribute 
% // it and/or modify it under the terms of the GNU General Public License 
% // as published by the Free Software Foundation, either version 3 of the 
% // License, or (at your option) any later version.
% //
% // Endmember Induction Algorithms Toolbox is distributed in the hope that it will
% // be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
% // of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% // General Public License for more details.
% //
% // You should have received a copy of the GNU General Public License
% // along with Endmember Induction Algorithms Toolbox. 
% // If not, see <http://www.gnu.org/licenses/>.
% // ====================================================================
%
%% [E,C] = EIA_NFINDR(data,p,maxit)
%
% Manuel Grana <manuel.grana[AT]ehu.es>
% Miguel Angel Veganzones <miguelangel.veganzones[AT]ehu.es>
% Grupo de Inteligencia Computacional (GIC), Universidad del Pais Vasco /
% Euskal Herriko Unibertsitatea (UPV/EHU)
% http://www.ehu.es/computationalintelligence
% 
% Copyright (2011) Grupo de Inteligencia Computacional @ Universidad del Pais Vasco, Spain.
% Copyright (2007) GRNPS group @ University of Extremadura, Spain. 
%
% N-FINDR endmembers induction algorithm.
% ------------------------------------------------------------------------------
% Input:   data      : column data matrix [nvariables x nsamples]
%          p         : number of endmembers to be induced. If not provided it is calculated by HFC method with tol=10^(-5)
%          maxit     : maximum number of iterations. Default = 3*p
%
% Output:  E         : set of induced endmembers [nvariables x p]
%          C         : induced endmembers indexes vector [nsamples] with {0,1} values, where '1' indicates that the corresponding sample has been identified as an endmember.
%
% Bibliographical references:
% [1] Winter, M. E., «N-FINDR: an algorithm for fast autonomous spectral end-member determination in hyperspectral data», presented at the Imaging Spectrometry V, Denver, CO, USA, 1999, vol. 3753, págs. 266-275.
function [E,C] = EIA_NFINDR(data,p,maxit)

%% Parameters
if (nargin < 1)
    error('Insufficient parameters');
end
if (nargin < 2 || isempty(p))
    p = EIA_HFC(data,10^(-5));
end
if (nargin < 3 || maxit <= 0)
    maxit = 3*p;
end

%% data size
[nvariables,nsamples] = size(data);

%% Dimensionality reduction by PCA
%[M,  ~, ~] = hyperMnf(data,7123,1);
[M,  ~, ~] = hyperMnf(data,191 ,28 );
%[M,  ~, ~] = hyperMnf(data,83 ,86 );
data_pca = M(1:p-1,:);

%% Initialization
E = zeros(nvariables,p);
C = zeros(1,nsamples);
IDX = zeros(1,p);
TestMatrix = zeros(p);
TestMatrix(1,:) = 1;
for i = 1:p
    idx = floor(rand*nsamples) + 1;
    TestMatrix(2:p,i) = data_pca(:,idx);
    IDX(i) = idx;
end
actualVolume = abs(det(TestMatrix)); % instead of: volumeactual = abs(det(MatrixTest))/(factorial(p-1));
it = 1;
v1 = -1;
v2 = actualVolume;

%% Algorithm
while it<=maxit && v2>v1
    for k=1:p
        for i=1:nsamples
            actualSample = TestMatrix(2:p,k);
            TestMatrix(2:p,k) = data_pca(:,i);
            volume = abs(det(TestMatrix));  % instead of: volume = abs(det(MatrixTest))/(factorial(p-1));
            if volume > actualVolume
                actualVolume = volume;
                IDX(k) = i;
            else
                TestMatrix(2:p,k) = actualSample;
            end
        end
    end
    it = it+1;
    v1 = v2;
    v2 = actualVolume;
end
for i = 1:p
    E(:,i) = data(:,IDX(i));
    C(IDX(i)) = 1;
end
