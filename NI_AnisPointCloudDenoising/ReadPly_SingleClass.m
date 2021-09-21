function [point_cloud] = ReadPly_SingleClass(ply_file)
% This function is for reading and extracting the point cloud matrix from 
% a PLY file
%
% INPUTS:  
%   ply_file   : a string containing the name of the PLY file to be opened.
%                This file contains the point cloud with only spatial
%                coordinate information.
% OUTPUTS:
%   point_cloud: N * 3 matrix in single precision.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (C) 2018-2020 Noiseless Imaging Oy - Tampere, Finland
%% Zhongwei Xu, Alessandro Foi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(ply_file);

% Extract Header
header = {};
header = [header;{fgetl(fid)}];
while (~strcmp(header{end},'end_header'))
    header = [header;{fgetl(fid)}];
end

% Extract Point Cloud Coordinates
point_cloud = single(cell2mat(textscan(fid,'%f %f %f')));
fclose(fid);
end