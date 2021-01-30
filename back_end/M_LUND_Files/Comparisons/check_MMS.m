function mms_in = check_MMS()

% Check if Partition Stability is in path:
if ispc 
  windows = 1;
else
  windows = 0;
end

pathCell = regexp(path, pathsep, 'split')';
n_dir = length(pathCell);

PartitionStability_in = 0;
GraphBasedClustering_in = 0;

for i = 1:n_dir
    
    full_path = pathCell{i};
    folder = full_path(max(strfind(full_path,'/'))+1:end);
    if windows
        PartitionStability_in = PartitionStability_in + strcmpi('PartitionStability', folder);
        GraphBasedClustering_in = GraphBasedClustering_in + strcmpi('GraphBasedClustering', folder);
    else
        PartitionStability_in = PartitionStability_in + strcmp('PartitionStability', folder);
        GraphBasedClustering_in = GraphBasedClustering_in + strcmp('GraphBasedClustering', folder);
    end
end

 mms_in = and(PartitionStability_in == 1, GraphBasedClustering_in == 1);

end