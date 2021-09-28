function [tot] = computeTot(h_star,count,LCS,tot,index,sigma_pcd)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

a=1:count;

maxh = max(abs(LCS(:,3)))/3;
for h_s = 1:length(h_star)
   
    for p=2:length(LCS)
        
            switch h_s
                case 1
                    if (LCS(p,1)<=h_star(h_s) && LCS(p,1)>=0 && LCS(p,2) <= h_star(h_s) &&  LCS(p,2)>=0 && abs(LCS(p,3)) <= max(6*sigma_pcd,h_star(h_s)/3)) 
                        if p<=index
                            tot(p-1) = tot(p-1)+1;
                        else
                            tot(p) = tot(p)+1;
                        end

                    end

                case 2
                    if (LCS(p,1)>=-h_star(h_s) && LCS(p,1)<=0 &&  LCS(p,2) <= h_star(h_s) &&  LCS(p,2)>=0 && abs(LCS(p,3)) <= max(6*sigma_pcd,h_star(h_s)/3))
                         if p<=index
                            tot(p-1) = tot(p-1)+1;
                        else
                            tot(p) = tot(p)+1;
                         end
                    end

                case 3
                    if (LCS(p,1)>=-h_star(h_s) && LCS(p,1)<=0 &&  LCS(p,2) >= -h_star(h_s) &&  LCS(p,2)<=0 && abs(LCS(p,3)) <= max(6*sigma_pcd,h_star(h_s)/3))
                        if p<=index
                            tot(p-1) = tot(p-1)+1;
                        else
                            tot(p) = tot(p)+1;
                        end
                    end

                case 4
                    if (LCS(p,1)<=h_star(h_s) && LCS(p,1)>=0 &&  LCS(p,2) >=-h_star(h_s) &&  LCS(p,2)<=0 && abs(LCS(p,3)) <= max(6*sigma_pcd,h_star(h_s)/3)) 
                         if p<=index
                            tot(p-1) = tot(p-1)+1;
                        else
                            tot(p) = tot(p)+1;
                        end
                    end
            end
    end
end
end



