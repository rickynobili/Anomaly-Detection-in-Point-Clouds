clear all;


pc = ReadPly_SingleClass('Bunny_noisy0_4delta.ply');
%pcd_ori = ReadPly_SingleClass('Bunny_ori.ply');

%pc = pointCloud(xyzPoints);

disp('Estimating the Noise Standard Deviation & Surface Sample Density ...');
tic; 
[sigma_pcd, dens_pcd] = pcd_stdEst_SingleClass(pc);
toc;   
fprintf('sigma_pcd =  %.5f;\n',sigma_pcd);
fprintf('dens_pcd =  %.5f;\n',dens_pcd);

pc = pcread('Bunny_noisy0_4delta.ply');

count = pc.Count;
stepsize = 1;
K = 800;

for index = 1000:stepsize:1000 %

    p = select(pc,index);
    
    KNN = findNearestNeighbors(pc,p.Location,K);
    
    points = zeros(numel(KNN),3);
    
    for i = 1:numel(KNN)
        p_i = select(pc,KNN(i));
        points(i,:) = p_i.Location;
    end
    
    npoints = zscore(points);
        
    coeff = pca(npoints);
        
    tcoeff = transpose(coeff);
        
    %biplot(coeff);
        
    LCS = zeros(count,3);
    
    LCS(1,:) = [ 0 0 0 ];
        
    for j = 1:stepsize:count
        
        if(j<index)
            pj = select(pc,j);

            temp = transpose(pj.Location - p.Location);

            LCS(j+1,:) = tcoeff*temp;
        end
        
        if(j>index)
            pj = select(pc,j);

            temp = transpose(pj.Location - p.Location);

            LCS(j,:) = tcoeff*temp;
        end
            
           
    end
        
    sharparam=-1;              % -1 zero order 0 first order (no sharpening) >0 sharpening
    gammaICI=1.05;             % ICI Gamma threshold
    directional_resolution=4;  % number of directions
    fusing=1;                  % fusing type   (1 classical fusing, 2 piecewise regular)
    addnoise=1;                % add noise to observation
    sigma_noise=0.1 ;
    
    h1=[5 10 15 20 25 30 ];
    %h2=max(1,ceil(h1*tan(0.5*pi/directional_resolution)));  % row vectors h1 and h2 need to have the same lenght
    %h2=ones(size(h1));
    h2=h1;
    lenh=length(h1);
    
    sig_winds=[ones(size(h1)); ones(size(h2))];    % Gaussian parameter
    beta=1;
    
    window_type=1;
    
    TYPE=11;            % TYPE IS A SYMMETRY OF THE WINDOW
    
    m = [0 0];
        
    %[kernels, kernels_higher_order] = function_CreateLPAKernels(m,h1,h2,TYPE,window_type,directional_resolution,sig_winds,beta);
        
    %scatter3(LCS(:,1),LCS(:,2),LCS(:,3));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    THETA = [pi/4 3*pi/4 5*pi/4 7*pi/4];
    
    order = 2;
    
    for quadrante=1:length(THETA)
    
        % trova indice rispetto a LCS dei punti nel directional neighb
        indice_z = 1; % serve per avere il vettore z con le varie stime di z per ogni valore di h considerato
        for h=h1
            pos = 1;
            for a=1:length(LCS)
                if (quadrante ==1) && (LCS(a,1)<h && LCS(a,1)>=0 && LCS(a,2) < h && LCS(a,2)>=0 && LCS(a,3) <= max(6*sigma_pcd,2*h)) %sto considerando il primo quadrante
                    indici(pos) = a;
                    pos = pos+1;
                end
                
                if (quadrante ==2) && (LCS(a,1)>-h && LCS(a,1)<=0 && LCS(a,2) < h && LCS(a,2)>=0 && LCS(a,3) <= max(6*sigma_pcd,2*h)) %sto considerando il primo quadrante
                    indici(pos) = a;
                    pos = pos+1;
                end
                
                if (quadrante ==3) && (LCS(a,1)>-h && LCS(a,1)<=0 && LCS(a,2) > -h && LCS(a,2)<=0 && LCS(a,3) <= max(6*sigma_pcd,2*h)) %sto considerando il primo quadrante
                    indici(pos) = a;
                    pos = pos+1;
                end
                
                if (quadrante ==4) && (LCS(a,1)<h && LCS(a,1)>=0 && LCS(a,2) >-h && LCS(a,2)<=0 && LCS(a,3) <= max(6*sigma_pcd,2*h)) %sto considerando il primo quadrante
                    indici(pos) = a;
                    pos = pos+1;
                end
            end

            PHI = [];
            
            for i1=0:order
                i2=0;
               while i1+i2 <= order
                   

                    for k=1:length(indici)
                        
                
                           phi(k) = (LCS(indici(k),1)^i1) * (LCS(indici(k),2)^i2);
                        
                    end

                    phi_t = transpose(phi);

                    if size(PHI)==0
                        PHI = phi_t;
                    else
                        PHI = [PHI,phi_t];
                    end

                    i2 = i2+1;

               end
            end

           for k=1:length(indici)
               g(k) = PHI(k,:)*pinv(transpose(PHI)*PHI)*transpose(PHI(1,:));
           end

           z(indice_z) = 0;

            for k=1:length(indici)
                z(indice_z) = z(indice_z) + LCS(k,3)*g(k);
            end
            indice_z = indice_z + 1;
        end

        l2z = norm(g,2);

        sigma_teta_h = sigma_pcd*(l2z^2);   

        for indice_z = 1:length(h1)
            D(indice_z,1) = z(indice_z)-gammaICI*sigma_teta_h;
            D(indice_z,2) = z(indice_z)+gammaICI*sigma_teta_h;
        end

        inf = D(1,1);
        sup = D(1,2);

        for indice_d=1:length(h1)-1

            ni = D(indice_d +1 , 1);
            ns= D(indice_d +1 , 2);

            if inf < ni && ni < sup < ns
                inf = ni;
                h_star(quadrante) = h1(indice_d);
            end

            if (sup > ns && ni < inf < ns)
                sup = ns;
                h_star(quadrante) = h1(indice_d);
            end

            if (inf > ns) || (sup < ni )
                h_star(quadrante) = h1(indice_d);
                break;
            end

            if (inf < ni) && (sup < ns )
                inf = ni;
                sup = ns;
                h_star(quadrante) = h1(indice_d);
            end

            if(indice_d == length(h1)-1)
                 h_star(quadrante) = h1(indice_d+1);
                 break;
            end

        end

    end
end
