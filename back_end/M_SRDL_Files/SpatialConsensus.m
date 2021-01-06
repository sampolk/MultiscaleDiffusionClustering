function [y,frac]=SpatialConsensus(Labels,Idx,M,N,I,J,WindowSize)

% For a set of labels, Labels, and a labelled point given by Idx, find the
% consensus label in a window of size WindowSize.

% Labels should be as a matrix

%Find spatial neighbors
Neighbors=FindNeighbors([I(Idx),J(Idx)],WindowSize,M,N);

%Neighbors=sub2ind([M,N],Neighbors(1,:),Neighbors(2,:));

%Determine labels of spatial neighbors

try
    %{
    for j=1:size(Neighbors,2)
        LocalLabels(j)=Labels(Neighbors(1,j),Neighbors(2,j));
    end
    %}
    LocalLabels=Labels(Neighbors(1,:),Neighbors(2,:));
    LocalLabels=LocalLabels(LocalLabels>0);
    
catch
    keyboard
end

if isempty(LocalLabels)
    y=0;
    frac=0;
else
    %y=mode(LocalLabels');
    y=mode(LocalLabels);
    frac=sum(LocalLabels==y)/length(LocalLabels);
end

end