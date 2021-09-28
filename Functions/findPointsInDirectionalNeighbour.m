function [indici] = findPointsInDirectionalNeighbour(quadrante,h,LCS,sigma_pcd)
indici=[];
pos = 1;
a=1:length(LCS);

maxh = max(abs(LCS(:,3)))/3;

for p=1 : length(LCS)
        
        switch quadrante
            
            case 1
                
                if(LCS(p,1)<h && LCS(p,1)>=0 && LCS(p,2) < h && LCS(p,2) >=0 && abs(LCS(p,3)) <= max(6*sigma_pcd,h/3)) %sto considerando il primo quadrante
                    indici(pos) = a(p);
                    pos = pos+1;
                end
                
            case 2

                if (LCS(p,1)>-h && LCS(p,1)<=0 && LCS(p,2) < h &&LCS(p,2)>=0 && abs(LCS(p,3)) <= max(6*sigma_pcd,h/3)) %sto considerando il primo quadrante
                    indici(pos) = a(p);
                    pos = pos+1;
                end
                
            case 3
                
                if (LCS(p,1)>-h && LCS(p,1)<=0 &&LCS(p,2)> -h && LCS(p,2)<=0 && abs(LCS(p,3)) <= max(6*sigma_pcd,h/3)) %sto considerando il primo quadrante
                    indici(pos) = a(p);
                    pos = pos+1;
                end
                
            case 4
                
                if (LCS(p,1)<h && LCS(p,1)>=0 && LCS(p,2) >-h && LCS(p,2)<=0 && abs(LCS(p,3)) <= max(6*sigma_pcd,h/3)) %sto considerando il primo quadrante
                    indici(pos) = a(p);
                    pos = pos+1;
                end
        end
end
end



