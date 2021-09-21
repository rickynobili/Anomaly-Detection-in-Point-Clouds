%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (C) 2018-2020 Noiseless Imaging Oy - Tampere, Finland
%% Zhongwei Xu, Alessandro Foi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = write_ply_only_pos(pcd, filename)
% This function is for exporting the point cloud matrix into a PLY file
% 
% INPUTS:  
%   pcd      : N * 3 matrix in single precision, contains the point cloud 
%              with only spatial coordinate information.
%   filename : a string of the name of the PLY file to be saved.
% OUTPUTS:
%   none 

fid = fopen(filename, 'w');

fprintf(fid, 'ply\n');
fprintf(fid, 'format ascii 1.0\n');
fprintf(fid, 'comment Noiseless Imaging Oy (Ltd) generated\n');
fprintf(fid, 'element vertex %d\n', size(pcd, 1));
fprintf(fid, 'property float x\n');
fprintf(fid, 'property float y\n');
fprintf(fid, 'property float z\n');
fprintf(fid, 'end_header\n');
fprintf(fid, '%f %f %f\n', pcd');

fclose(fid);

end