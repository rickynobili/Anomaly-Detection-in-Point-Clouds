%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (C) 2018-2020 Noiseless Imaging Oy - Tampere, Finland
%% Zhongwei Xu, Alessandro Foi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo program to implement the algorithm proposed in: 
% Z. Xu and A. Foi, "Anisotropic Denoising of 3D Point Clouds by 
% Aggregation of Multiple Surface-Adaptive Estimates", 
% IEEE Transactions on Visualization and Computer Graphics, 2019.
% http://doi.org/10.1109/TVCG.2019.2959761

%% add folder and subfolders to path 

addpath(genpath('.'));

%% Read Noisy and Ground-Truth Point Clouds in single precision
%  more testing noisy point clouds and its corresponding ground truth are 
%  in folders:
%    './Noisy_point_clouds'  
%    './NoiseFree_point_clouds'.

% pcd_noi = ReadPly_SingleClass('Dodecahedron_noi.ply');
% pcd_ori = ReadPly_SingleClass('Dodecahedron_ori.ply');

% pcd_noi = ReadPly_SingleClass('cube49_noisy0_4delta.ply');
% pcd_ori = ReadPly_SingleClass('Cube49_ori.ply');

% pcd_noi = ReadPly_SingleClass('Sphere_noisy0_4delta.ply');
% pcd_ori = ReadPly_SingleClass('Sphere_ori.ply');

% pcd_noi = ReadPly_SingleClass('Fandisk_noisy0_8delta.ply');
% pcd_ori = ReadPly_SingleClass('Fandisk_ori.ply');

pcd_noi = ReadPly_SingleClass('Bunny_noisy0_8delta.ply');
pcd_ori = ReadPly_SingleClass('Bunny_ori.ply');

%pcd_noi = ReadPly_SingleClass('Armadillo_noisy0_8delta.ply');
%pcd_ori = ReadPly_SingleClass('Armadillo_ori.ply');

% pcd_noi = ReadPly_SingleClass('Shutter_blind.ply');

% pcd_noi = ReadPly_SingleClass('Iron.ply');

%% Estimate Noise Standard Deviation & Surface Sample Density (not compiled)
disp('Estimating the Noise Standard Deviation & Surface Sample Density ...');
tic; 
[sigma_pcd, dens_pcd] = pcd_stdEst_SingleClass(pcd_noi);
toc;   
fprintf('sigma_pcd =  %.5f;\n',sigma_pcd);
fprintf('dens_pcd =  %.5f;\n',dens_pcd);

%% Denoising
disp('Denoising ...');
tstart = tic; 
[pcd_de] = PCD_Filtering_LPAICI_mex(pcd_noi,sigma_pcd,dens_pcd);

t1 = toc(tstart);
disp(['Run time is ',num2str(t1,'%.02f'),' seconds']);

%%
figure; 
subplot(1,2,1); plot3(pcd_noi(:,1),pcd_noi(:,2),pcd_noi(:,3),'.'); 
axis equal;
title('Noisy');
subplot(1,2,2); plot3(pcd_de(:,1),pcd_de(:,2),pcd_de(:,3),'.'); 
axis equal;
title('Denoised');

%% compute Mean squared point-to-surface distance (not compiled)
RSMSE_p2s = SquareRoot_MeanPoint2SurfError_SingleClass(pcd_de,pcd_ori); 
fprintf('Root mean squared point-to-surface distance = %.5f;\n',RSMSE_p2s);

%% Save denoised point cloud
write_ply_only_pos(pcd_de, 'pcd_denoised.ply');
fprintf('Denoised point cloud has been saved in current file directory\n');

