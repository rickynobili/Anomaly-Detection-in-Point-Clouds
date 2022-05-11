clear all;
%% Read the .ply file
%x = load("reducedPotato021.mat");
%pc=x.reducedpotato021;
%pc = pointCloud(pc);
%pc = x.pc;
%pc = ReadPly_SingleClass('lastra_ondulata_2_skretches_50x50_zup.ply'); %Bunny_noisy0_2delta.ply %Fandisk_noisy0_8delta.ply %Fandisk_noisy0_4delta.ply
pc = pcread('damaged_column1_2.ply');
loc = single(pc.Location);
rng(2022);
%pcd_ori = ReadPly_SingleClass('Bunny_ori.ply');

%% Add Gaussian noise to the point cloud
var =0.001;
m = 0;
noise = normrnd(m,var,length(loc),3);
loc = loc + noise;
pc = pointCloud(loc,'Color',pc.Color);

%pc = pcdenoise(pc);
%num_points = 60000;
%pc = pcdownsample(pc,'random', num_points/length(pc.Location));
loc = single(pc.Location);

%% Estimation of the Noise Standard Deviation & Surface Sample Density 

disp('Estimating the Noise Standard Deviation & Surface Sample Density ...');
tic; 
[sigma_pcd, dens_pcd] = pcd_stdEst_SingleClass(loc);
toc;   
fprintf('sigma_pcd =  %.5f;\n',sigma_pcd);
fprintf('dens_pcd =  %.5f;\n',dens_pcd);


clear location noise;

%% Initialize useful variables


count = pc.Count;               % total number of points in the point cloud
stepsize = 3;                       % size of the step when considering the points in the point cloud

global K;
K=80;                                  % KNN

MOD=2;                              % useful to print which points have been analysed

tot = zeros(1,count);          % initialization of the variable that specifies the number of directional neighbours each point falls into

THETA = [pi/4 3*pi/4 5*pi/4 7*pi/4];

order = 2;                           % order for the local polynomial model

gammaICI=2;                     % ICI Gamma threshold
directional_resolution=4;    % number of quadrants that have to be considered (solitamente tutti e 4 quelli identificati da assi x e y dell'LCS)

h1=[(3*(sqrt(2))^0)/sqrt(dens_pcd) (3*(sqrt(2))^1)/sqrt(dens_pcd) (3*(sqrt(2))^2)/sqrt(dens_pcd) (3*(sqrt(2))^3)/sqrt(dens_pcd) (3*(sqrt(2))^4)/sqrt(dens_pcd) (3*(sqrt(2))^5)/sqrt(dens_pcd) (3*(sqrt(2))^8)/sqrt(dens_pcd)];
%h1 = [0.1 0.2 0.3 0.5 1 2];

%% Define the real defects as the colored points in the point cloud

% pc.Color(pc.Color==254)=255;
% real_defects = zeros(1,count);
% 
% if not(isempty(pc.Color))
%     for p = 1:count
%         if (not(isequal(pc.Color(p,:),[255 255 255])))
%             real_defects(p) = 1;
%         end
%     end
% end


%% Execution of the algorithm

for index = 1:stepsize:count  % consider each point in the point cloud
    
    if mod(index,MOD)==0
        X = ['index ', num2str(index), ' analized'];
        disp(X);
    end
    
    clear h_star indici KNN LCS g PHI points;
    
    h_star = zeros(1,4);    % initialize h_star

    p = select(pc,index);   % select a point in the point cloud
    
    LCS = CreateLCS(p,pc,index,K);  % create the LCS
    
    for quadrante=1:length(THETA)   % consider each quadrant
    
        % trova indice rispetto a LCS dei punti nel directional neighb
        indice_z = 1; % serve per avere il vettore z con le varie stime di z per ogni valore di h considerato
        
        for h=h1
            
            % find the indices wrt the LCS of the points 
            % that fall in the considered directional neighbour
            indici = findPointsInDirectionalNeighbour(quadrante,h,LCS,sigma_pcd);   
            
            
            % compute the matrix PHI
            PHI = computePHI(order,indici,LCS);
            
            % Compute g
            
            
            k=1:length(indici);
            g = PHI(k,:)*pinv(transpose(PHI)*PHI)*transpose(PHI(1,:));
            
          
           % Compute z
           z(indice_z) = 0;

            for k=1:length(indici)
                z(indice_z) = z(indice_z) + LCS(indici(k),3)*g(k);
            end
            indice_z = indice_z + 1;
        end
        
        % ICI
        l2g = norm(g,2)^2;

        sigma_teta_h = (sigma_pcd^2)*l2g;
        sigma_teta_h = sqrt(sigma_teta_h);
        
        h_star = performICI(h_star,sigma_teta_h,gammaICI,h1,quadrante,z);
        
    end
    
    % compute the number of DNs each point falls into
    tot = computeTot(h_star,count,LCS,tot,index,sigma_pcd);
    
end

%% Compute the metrics
% visualize the ROC curve and compute the AUC
anomalyScore = 1./tot;
anomalyScore(anomalyScore==inf) = 1;
[f,xi] = ksdensity(anomalyScore,anomalyScore);
toc;
%if(isempty(pc.Color))
    
%pcshow(pc.Location,transpose(anomalyScore));
anomaly_test = anomalyScore;
%scatter(xi,f,1);
[x,y,THR,AUC,OPT] = perfcurve(real_defects,anomalyScore,1);
ThresholdForOptROCpt = THR((x==OPT(1))&(y==OPT(2)))
anomaly_test(anomaly_test<(min(anomaly_test)+std(anomaly_test)))=0;
%anomaly_test(anomaly_test<ThresholdForOptROCpt)=0;

pcshow(pc.Location,transpose(anomaly_test));

found_defects = anomaly_test;
found_defects(found_defects > 0) = 1;

[precision,recall,accuracy,F1,FNR] = ComputeMetrics(real_defects,found_defects);
 
    
[x,y,THR,AUC,OPT] = perfcurve(real_defects,anomalyScore,1);
subplot(2,2,3);
plot(x,y,'LineWidth',4);
% % 
cmatrix = (transpose(found_defects) * [1 0 0]);
tmatrix = (transpose(real_defects) * [1 0 0]);

ptc = pointCloud(pc.Location,'Color',cmatrix);
ptt = pointCloud(pc.Location,'Color',tmatrix);

subplot(2,2,1)
pcshow(ptt);
subplot(2,2,2)
pcshow(ptc);

toc;
