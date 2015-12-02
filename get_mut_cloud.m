function [ cloud ] = get_mut_cloud( cloud )
% remove from the cloud all the positions that get onto the consensus, 
% leaving the mutations cloud
% input:
% clud = the original cloud
% output:
% cloud = the mutation cloud

[ ~, max_poses] = max(cloud);
for i=1:size(cloud,2)
    cloud(max_poses(i), i) = 0;
end

end

