function [indici] = findPointsInDirectionalNeighbour(quadrante,h,LCS,sigma_pcd)
indici=[];
pos = 1;
    for a=1:length(LCS)
        switch quadrante
            
            case 1
                
                if(LCS(a,1)<h && LCS(a,1)>=0 && LCS(a,2) < h && LCS(a,2)>=0 && abs(LCS(a,3)) <= max(6*sigma_pcd,2*h))/2 %sto considerando il primo quadrante
                    indici(pos) = a;
                    pos = pos+1;
                end
                
            case 2

                if (LCS(a,1)>-h && LCS(a,1)<=0 && LCS(a,2) < h && LCS(a,2)>=0 && abs(LCS(a,3)) <= max(6*sigma_pcd,2*h))/2 %sto considerando il primo quadrante
                    indici(pos) = a;
                    pos = pos+1;
                end
                
            case 3
                
                if (LCS(a,1)>-h && LCS(a,1)<=0 && LCS(a,2) > -h && LCS(a,2)<=0 && abs(LCS(a,3)) <= max(6*sigma_pcd,2*h))/2 %sto considerando il primo quadrante
                    indici(pos) = a;
                    pos = pos+1;
                end
                
            case 4
                
                if (LCS(a,1)<h && LCS(a,1)>=0 && LCS(a,2) >-h && LCS(a,2)<=0 && abs(LCS(a,3)) <= max(6*sigma_pcd,2*h))/2 %sto considerando il primo quadrante
                    indici(pos) = a;
                    pos = pos+1;
                end
        end
    end
end

