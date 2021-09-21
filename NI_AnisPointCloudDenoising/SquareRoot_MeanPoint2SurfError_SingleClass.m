function [SMSE_point_to_surface,distance_point_to_surface] = SquareRoot_MeanPoint2SurfError_SingleClass(point_cloud_check,point_cloud_ori)
% This program compares two point clouds (estimated one & ground truth), 
% which can have different number of points, by computing the sqaured mean
% point-to-surface distance between them.
%
% INPUTS:
%   POINT_CLOUD_CHECK: M * 3 matrix of the estimated point cloud
%   POINT_CLOUD_ORI  : N * 3 matrix of the ground truth point cloud
% OUTPUTS:
%   SMSE_POINT_TO_SURFACE    : sqaured mean point-to-surface distance
%   DISTANCE_POINT_TO_SURFACE: M * 1 matrix of point-to-surface distance
%                              for each point in POINT_CLOUD_CHECK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (C) 2018-2020 Noiseless Imaging Oy - Tampere, Finland
%% Zhongwei Xu, Alessandro Foi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_point_est = single(size(point_cloud_check,1));

distance_point_to_surface = single(zeros(size([1:num_point_est],2),1));

Idx = knnsearch(point_cloud_ori,point_cloud_check,'K',1,'NSMethod','kdtree');
knn_neigh = knnsearch(point_cloud_ori,point_cloud_ori,'K',5,'NSMethod','kdtree');

kk = 0;
for i = 1:num_point_est
	kk = kk+1;
    I = Idx(i);
    knn_neigh_I = point_cloud_ori(knn_neigh(I,:),:);
	[V_rotate_c,I_rotate_c] = local_coordinate_for_ICI(knn_neigh_I); 
	v_normal = V_rotate_c(:,I_rotate_c(3));
	
	distance_point_to_surface(kk) = norm((point_cloud_check(i,:) - point_cloud_ori(I,:))*v_normal);
end

SMSE_point_to_surface = sqrt(mean(distance_point_to_surface.^2));
end

function [V2, I, D1] = local_coordinate_for_ICI(knn)
%   LOCAL_COORDINATE function reads the KNN matrix of a point, 
%   and compute the principle components of it.

M = bsxfun(@minus,knn,mean(knn));
P1 = M'*M;
[V,D] = eig(P1,'vector');
V2 = real(V);
D1 = D(1);
[~,I] = sort(D,'descend');
end