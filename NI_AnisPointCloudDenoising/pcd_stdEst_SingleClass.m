function [sigma_pcd, dens_pcd] = pcd_stdEst_SingleClass(pcd)
% This functions estimates the standard deviation of isotropic additive 
% Gaussian noise in a point cloud, as well as the sampling density 
% of the point cloud.
%
% see
% Z. Xu and A. Foi, "Anisotropic Denoising of 3D Point Clouds by 
% Aggregation of Multiple Surface-Adaptive Estimates", 
% IEEE Transactions on Visualization and Computer Graphics, 2019.
% http://doi.org/10.1109/TVCG.2019.2959761
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (C) 2018-2020 Noiseless Imaging Oy - Tampere, Finland
% Zhongwei Xu, Alessandro Foi
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_k_stdEst = single(50); % define the K numbers of neighbors used for std estimate

num_p = single(size(pcd,1));
d3 = zeros(num_p,1,'single');
d12 = zeros(num_p,1,'single');

Idx = int32(knnsearch(pcd,pcd,'K',num_k_stdEst,'NSMethod','kdtree'));

for i = 1:num_p
    
    I = Idx(i,:);
    knn = pcd(I,:);
    
    [V, I_local] = local_coordinate_for_ICI(knn);   % computing Principle Component Vector
    knn_new = bsxfun(@minus, knn, knn(1,:)) * V;    % convert to LCS
    pcd_patch = [knn_new(:,I_local(1)),knn_new(:,I_local(2)),knn_new(:,I_local(3))];
    
    distance = sum(bsxfun(@minus,pcd_patch(:,1:2),pcd_patch(1,1:2)).^2,2);
    distance(1) = inf;
    
    [~,I2] = min(distance);
    z_nearest = pcd_patch(I2,3);
    d3(i) = abs(z_nearest - pcd_patch(1,3))/sqrt(2);
    d12(i) = var(pcd_patch(:,1:2)*[1 1i]')/size(pcd_patch,1);
end

sigma_pcd = median(d3)/0.6745;
dens_pcd = 1/(median(d12)*2*pi);

%============= recompute with bigger K, if sigma_pcd is large =============
if  (sigma_pcd * sqrt(dens_pcd) > 1.5 && sigma_pcd * sqrt(dens_pcd) < 3.5)
    num_k_stdEst = single(200); % define the K numbers of neighbors used for std estimate
    
    num_p = single(size(pcd,1));
    d3 = zeros(num_p,1,'single');
    d12 = zeros(num_p,1,'single');
    
    Idx = int32(knnsearch(pcd,pcd,'K',num_k_stdEst,'NSMethod','kdtree'));
    
    for i = 1:num_p
        I = Idx(i,:);
        knn = pcd(I,:);
        [V, I_local] = local_coordinate_for_ICI(knn);   % computing Principle Component Vector
        knn_new = bsxfun(@minus, knn, knn(1,:)) * V;    % convert to LCS
        pcd_patch = [knn_new(:,I_local(1)),knn_new(:,I_local(2)),knn_new(:,I_local(3))];
        
        distance = sum(bsxfun(@minus,pcd_patch(:,1:2),pcd_patch(1,1:2)).^2,2);
        distance(1) = inf;
        
        [~,I2] = min(distance);
        z_nearest = pcd_patch(I2,3);
        d3(i) = abs(z_nearest - pcd_patch(1,3))/sqrt(2);
        d12(i) = var(pcd_patch(:,1:2)*[1 1i]')/size(pcd_patch,1);
    end
    
    sigma_pcd = median(d3)/0.6745;
    dens_pcd = 1/(median(d12)*2*pi);
end

if sigma_pcd * sqrt(dens_pcd) >= 3.5 && sigma_pcd * sqrt(dens_pcd) < 4.5
    num_k_stdEst = single(300); % define the K numbers of neighbors used for std estimate
    
    num_p = single(size(pcd,1));
    d3 = zeros(num_p,1,'single');
    d12 = zeros(num_p,1,'single');
    
    Idx = int32(knnsearch(pcd,pcd,'K',num_k_stdEst,'NSMethod','kdtree'));
    
    for i = 1:num_p
        I = Idx(i,:);
        knn = pcd(I,:);
        [V, I_local] = local_coordinate_for_ICI(knn);   % computing Principle Component Vector
        knn_new = bsxfun(@minus, knn, knn(1,:)) * V;    % convert to LCS
        pcd_patch = [knn_new(:,I_local(1)),knn_new(:,I_local(2)),knn_new(:,I_local(3))];
        
        distance = sum(bsxfun(@minus,pcd_patch(:,1:2),pcd_patch(1,1:2)).^2,2);
        distance(1) = inf;
        
        [~,I2] = min(distance);
        z_nearest = pcd_patch(I2,3);
        d3(i) = abs(z_nearest - pcd_patch(1,3))/sqrt(2);
        d12(i) = var(pcd_patch(:,1:2)*[1 1i]')/size(pcd_patch,1);
    end
    
    sigma_pcd = median(d3)/0.6745;
    dens_pcd = 1/(median(d12)*2*pi);
end

if sigma_pcd * sqrt(dens_pcd) >= 4.5
    num_k_stdEst = single(500); % define the K numbers of neighbors used for std estimate
    
    num_p = single(size(pcd,1));
    d3 = zeros(num_p,1,'single');
    d12 = zeros(num_p,1,'single');
    
    Idx = int32(knnsearch(pcd,pcd,'K',num_k_stdEst,'NSMethod','kdtree'));
    
    for i = 1:num_p
        I = Idx(i,:);
        knn = pcd(I,:);
        [V, I_local] = local_coordinate_for_ICI(knn);   % computing Principle Component Vector
        knn_new = bsxfun(@minus, knn, knn(1,:)) * V;    % convert to LCS
        pcd_patch = [knn_new(:,I_local(1)),knn_new(:,I_local(2)),knn_new(:,I_local(3))];
        
        distance = sum(bsxfun(@minus,pcd_patch(:,1:2),pcd_patch(1,1:2)).^2,2);
        distance(1) = inf;
        
        [~,I2] = min(distance);
        z_nearest = pcd_patch(I2,3);
        d3(i) = abs(z_nearest - pcd_patch(1,3))/sqrt(2);
        d12(i) = var(pcd_patch(:,1:2)*[1 1i]')/size(pcd_patch,1);
    end
    
    sigma_pcd = median(d3)/0.6745;
    dens_pcd = 1/(median(d12)*2*pi);
end
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