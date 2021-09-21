clear all;
%% Read the .ply file

%pc = ReadPly_SingleClass('lastra_ondulata_2_skretches_50x50_zup.ply'); %Bunny_noisy0_2delta.ply %Fandisk_noisy0_8delta.ply %Fandisk_noisy0_4delta.ply
pc = pcread('lastra_rugosa_3_difetti_sconnessi_colored.ply');
loc = single(pc.Location);
%pcd_ori = ReadPly_SingleClass('Bunny_ori.ply');

%% Add Gaussian noise to the point cloud
var =0.01;
m = 0;
noise = normrnd(m,var,length(loc),3);
loc = loc + noise;
pc = pointCloud(loc,'Color',pc.Color);

pc = pcdenoise(pc);
%pc = pcdownsample(pc,'random', 9000/length(pc.Location));
loc = single(pc.Location);

%% Estimation of the Noise Standard Deviation & Surface Sample Density 

disp('Estimating the Noise Standard Deviation & Surface Sample Density ...');
tic; 
[sigma_pcd, dens_pcd] = pcd_stdEst_SingleClass(loc);
toc;   
fprintf('sigma_pcd =  %.5f;\n',sigma_pcd);
fprintf('dens_pcd =  %.5f;\n',dens_pcd);

% Recreate the pointcloud with the added noise and downsampling


clear location noise;

%% Initialize useful variables


count = pc.Count;               % total number of points in the point cloud
stepsize = 1;                       % size of the step when considering the points in the point cloud

global K;
K=80;                                  % KNN

MOD=2;                              % useful to print which points have been analysed

tot = zeros(1,count);          % initialization of the variable that specifies the number of directional neighbours each point falls into

THETA = [pi/4 3*pi/4 5*pi/4 7*pi/4];

order = 2;                           % order for the local polynomial model

gammaICI=1.15;                     % ICI Gamma threshold

directional_resolution=4;    % number of directions

h1=[(3*(sqrt(2))^0)/sqrt(dens_pcd) (3*(sqrt(2))^1)/sqrt(dens_pcd) (3*(sqrt(2))^2)/sqrt(dens_pcd) (3*(sqrt(2))^3)/sqrt(dens_pcd) (3*(sqrt(2))^4)/sqrt(dens_pcd) (3*(sqrt(2))^5)/sqrt(dens_pcd)];
%h1 = [0.2 0.5 0.8 1]

%% Define the real defects as the colored points in the point cloud

real_defects = zeros(1,count);

for p = 1:count
	if (not(isequal(pc.Color(p,:),[255 255 255])))
		real_defects(p) = 1;
	end
end

%% Execution of the algorithm

for index = 1:stepsize:count  % consider each point in the point cloud
    
    if mod(index,MOD)==0
        X = ['index ', num2str(index), ' analized'];
        disp(X);
    end
    
    clear h_star indici KNN LCS g PHI points;
    
    h_star = zeros(1,4);    % initialize h_star

    p = select(pc,index);   % select a point in the point cloud
    
    [LCS,tcoef] = CreateLCS(p,pc,index,K);  % create the LCS
    
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
        
%         D = zeros(length(h1),2);
%         
%         % Compute the intervals for the ICI
%         for indice_z = 1:length(h1)
%             D(indice_z,1) = z(indice_z)-gammaICI*sigma_teta_h;
%             D(indice_z,2) = z(indice_z)+gammaICI*sigma_teta_h;
%         end
% 
%         inferiore = D(1,1);
%         superiore = D(1,2);
%         
%         for indice_d=1:length(h1)-1
% 
%             ni = D(indice_d +1 , 1);
%             ns= D(indice_d +1 , 2);
% 
%             if (inferiore < ni && ni < superiore && superiore < ns)
%                 inferiore = ni;
%                 h_star(quadrante) = h1(indice_d+1);
%             end
% 
%             if (superiore > ns && ni < inferiore && inferiore < ns)
%                 superiore = ns;
%                 h_star(quadrante) = h1(indice_d+1);
%             end
% 
%             if (inferiore > ns) || (superiore < ni )
%                 h_star(quadrante) = h1(indice_d);
%                 break;
%             end
% 
%             if (inferiore < ni) && (superiore > ns )
%                 inferiore = ni;
%                 superiore = ns;
%                 h_star(quadrante) = h1(indice_d+1);
%             end
% 
%             if(indice_d == length(h1)-1)
%                  h_star(quadrante) = h1(indice_d+1);
%                  break;
%             end
%         end
    end
    
    tot = computeTot(h_star,count,LCS,tot,index,sigma_pcd);
    
%         for h_s = 1:length(h_star)
%             for a=1:count
%                 if (h_s==1) && (LCS(a,1)<=h_star(h_s) && LCS(a,1)>=0 && LCS(a,2) <= h_star(h_s) && LCS(a,2)>=0 && abs(LCS(a,3)) <= max(3*sigma_pcd,h_star(h_s)))
%                     if a==1
%                       tot(index) = tot(index) +1;
% 
%                     elseif a<= index
%                         tot(a-1) = tot(a-1)+1;
%                     elseif a>index
%                         tot(a) = tot(a)+1;
%                     end
%                 end
% 
%                 if (h_s ==2) && (LCS(a,1)>=-h_star(h_s) && LCS(a,1)<=0 && LCS(a,2) <= h_star(h_s) && LCS(a,2)>=0 && abs(LCS(a,3)) <= max(3*sigma_pcd,h_star(h_s)))
%                     if a==1
%                       tot(index) = tot(index) +1;
% 
%                     elseif a<= index
%                         tot(a-1) = tot(a-1)+1;
% 
%                     elseif a>index
%                         tot(a) = tot(a)+1;
%                     end
%                 end
% 
%                 if (h_s ==3) && (LCS(a,1)>=-h_star(h_s) && LCS(a,1)<=0 && LCS(a,2) >= -h_star(h_s) && LCS(a,2)<=0 && abs(LCS(a,3)) <= max(3*sigma_pcd,h_star(h_s)))
%                     if a==1
%                       tot(index) = tot(index) +1;
% 
%                     elseif a<= index
%                         tot(a-1) = tot(a-1)+1;
% 
%                     elseif a>index
%                         tot(a) = tot(a)+1;
%                     end
%                 end
% 
%                 if (h_s ==4) && (LCS(a,1)<=h_star(h_s) && LCS(a,1)>=0 && LCS(a,2) >=-h_star(h_s) && LCS(a,2)<=0 && abs(LCS(a,3)) <= max(3*sigma_pcd,h_star(h_s))) 
%                     if a==1
%                       tot(index) = tot(index) +1;
% 
%                     elseif a<= index
%                         tot(a-1) = tot(a-1)+1;
% 
%                     elseif a>index
%                         tot(a) = tot(a)+1;
%                     end       
%                 end
%             end
%         end
end

% visualize the ROC curve and compute the AUC
anomalyScore = 1./tot;

[x,y,THR,AUC,OPT] = perfcurve(real_defects,anomalyScore,1);
subplot(2,2,3);
plot(x,y,'LineWidth',4);


cmatrix = (transpose(tot) * [1 0 0]);
tmatrix = (transpose(real_defects) * [0 1 0]);

% Cerca il valore start tale che se un punto rientra in meno di start
% directional neighbours allora è considerato anomalia e cerca start
% massimizzando la F1

start = round(mean(tot));

found_defects = findDefects(start,cmatrix,tmatrix,pc);

[precision,recall,accuracy,F1,FPR] = ComputeMetrics(real_defects,found_defects);

F1t = F1;

while(F1t>=F1)
    F1 = F1t;
    start = start - 1;
    found_defects = findDefects(start,cmatrix,tmatrix,pc);
    [precision,recall,accuracy,F1t,FPR] = ComputeMetrics(real_defects,found_defects);
end

if(start == round(mean(tot))-1)
    
    start = round(mean(tot));

    found_defects = findDefects(start,cmatrix,tmatrix,pc);

    [precision,recall,accuracy,F1,FPR] = ComputeMetrics(real_defects,found_defects);
    
    F1t = F1;

    while(F1t>=F1)
        F1 = F1t;
        start = start + 1;
        found_defects = findDefects(start,cmatrix,tmatrix,pc);
        [precision,recall,accuracy,F1t,FPR] = ComputeMetrics(real_defects,found_defects);
    end
end
toc;