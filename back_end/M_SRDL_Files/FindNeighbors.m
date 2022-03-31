%For a pixel coordinates in a matrix, find the size^2-1
%spatial nearest neighbors.

function Neighbors=FindNeighbors(Coordinates, Size, m, n)

Neighbors=[];

for i=1:2*Size+1
    for j=1:2*Size+1
        if ( (Coordinates(1)-Size-1+i)>0 && (Coordinates(2)-Size-1+j)>0 && (Coordinates(1)-Size-1+i)<(m+1) && (Coordinates(2)-Size-1+j)<(n+1))
            Neighbors=[Neighbors,[Coordinates(1)-Size-1+i,Coordinates(2)-Size-1+j]'];
        end
    end
end


